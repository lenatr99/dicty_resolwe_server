""".. Ignore pydocstyle D400.

===========================
Expose User & CSRF on API
===========================

"""

import logging

import uuid
from django.utils import timezone
from rest_framework.decorators import action

import django_filters as filters
from django.utils.decorators import method_decorator
from django.middleware.csrf import get_token
from django.views.decorators.csrf import ensure_csrf_cookie

from rest_framework import exceptions, mixins, viewsets
from rest_framework.permissions import AllowAny, IsAdminUser, IsAuthenticated
from resolwe_bio.kb.views import FeatureViewSet
from rest_framework.response import Response
from resolwe.flow.views import DataViewSet as BaseDataViewSet  # import from the package
from resolwe_bio.kb.models import Feature

from django.contrib.auth import get_user_model

from resolwe.flow.filters import OrderingFilter
from resolwe.flow.models import Entity as Sample
# views.py
from rest_framework.viewsets import ViewSet
from django.db.models.functions import Coalesce
from django.db.models import Value, Q

from resolwe.flow.models import Data
from django.core.cache import cache

from .serializers import UserSerializer, BasketCreateSerializer

logger = logging.getLogger(__name__)
User = get_user_model()


class UserFilter(filters.FilterSet):
    username = filters.CharFilter(field_name="username", lookup_expr="icontains")
    email = filters.CharFilter(field_name="email", lookup_expr="icontains")
    is_staff = filters.BooleanFilter(field_name="is_staff")

    class Meta:
        model = User
        fields = ["username", "email", "is_staff"]


class UserViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    """User endpoint.

    - `GET /api/user/?current_only=1` → returns the currently authenticated user.
    - `GET /api/user/` → staff only (lists users, filterable/orderable).
    """

    queryset = User.objects.all().order_by("id")
    serializer_class = UserSerializer
    filter_backends = [filters.rest_framework.DjangoFilterBackend, OrderingFilter]
    filterset_class = UserFilter
    ordering_fields = ("id", "username", "email", "is_staff")

    def list(self, request, *args, **kwargs):
        if request.query_params.get("current_only") == "1":
            # Any authenticated user may fetch themselves.
            serializer = self.get_serializer(request.user)
            return Response(serializer.data)

        # Otherwise, only staff may list users.
        if not request.user.is_staff:
            raise exceptions.PermissionDenied("Listing users requires staff privileges.")

        return super().list(request, *args, **kwargs)


class CsrfTokenViewSet(viewsets.ViewSet):
    """CSRF token endpoint.

    - `GET /api/base/csrftoken` → returns `{ "csrfToken": "<token>" }`
      and ensures the `csrftoken` cookie is set.
    """

    permission_classes = (AllowAny,)

    @method_decorator(ensure_csrf_cookie)
    def list(self, request):
        return Response({"csrfToken": get_token(request)})



class DifferentialExpressionViewSet(BaseDataViewSet):
    http_method_names = ["get"]  # read-only

    def list(self, request):
        response = {}
        return Response(response)

class BasketViewSet(mixins.CreateModelMixin, viewsets.GenericViewSet):
    serializer_class = BasketCreateSerializer
    http_method_names = ["post", "options", "head"]

    @action(detail=False, methods=["post"], url_path=r"_/add_samples")
    def add_samples(self, request):
        # reuse your create() logic in here (or call self.create(request))
        ser = self.get_serializer(data=request.data)
        ser.is_valid(raise_exception=True)
        v = ser.validated_data

        found = Sample.objects.filter(id__in=v["samples"])
        found_ids = set(found.values_list("id", flat=True))
        ignored = [i for i in v["samples"] if i not in found_ids]

        rows = (
            Data.objects
            .filter(process__type__startswith="data:expression", entity__id__in=found_ids)
            .values_list("output__species", "input__species",
                        "output__source",  "input__source", "input__gene_id_database")
        )

        orgs = {a for a, b, *_ in rows if a} | {b for a, b, *_ in rows if b}
        srcs = {s for *_, s1, s2, s3 in rows for s in (s1, s2, s3) if s}


        basket_id = str(uuid.uuid4())
        selected_ids = [s.id for s in found]  # final sample IDs
        cache.set(f"basket:{basket_id}", selected_ids, 86400)  # 1 day TTL

        return Response({
            "id": basket_id,
            "modified": timezone.now().isoformat(),
            "ignored": ignored,
            "duplicated": [],
            "permitted_organisms": list(orgs),
            "permitted_sources": list(srcs),
            "conflict_organisms": [],
            "conflict_sources": [],
        })
    
class BasketExpressionsViewSet(viewsets.ViewSet):
    def list(self, request):
        basket = request.query_params.get("basket", "")
        ids = cache.get(f"basket:{basket}", [])
        if not ids:
            return Response([])

        rows = (
            Data.objects
            .filter(process__type__startswith="data:expression", entity_id__in=ids)
            .order_by("id")                                # <- ensure ascending order
            .values("id", "input", "output")
        )

        return Response([
            {
                "id": r["id"],
                "exp_type": (r["output"] or {}).get("exp_type")
                           or (r["input"] or {}).get("exp_type")
                           or "polyA"
            }
            for r in rows
        ])






class CustomFeatureViewSet(FeatureViewSet):
    """Custom FeatureViewSet that returns only results without pagination wrapper."""
    
    def list(self, request, *args, **kwargs):
        queryset = self.filter_queryset(self.get_queryset())
        serializer = self.get_serializer(queryset, many=True)
        return Response(serializer.data)


class GeneListByIdsViewSet(FeatureViewSet):
    def create(self, request, *args, **kwargs):
        ids = {s.strip() for s in request.data.get("gene_ids", []) if s.strip()}
        qs = self.get_queryset().filter(Q(feature_id__in=ids) | Q(name__in=ids) | Q(aliases__overlap=list(ids)))
        src = request.data.get("source"); sp = request.data.get("species")
        if src: qs = qs.filter(source__iexact=src)
        if sp:  qs = qs.filter(species__iexact=sp)
        page = self.paginate_queryset(qs); ser = self.get_serializer(page or qs, many=True)
        return self.get_paginated_response(ser.data) if page is not None else Response(ser.data)

