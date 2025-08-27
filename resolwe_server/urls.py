from django.urls import include, path
from django.contrib.auth import authenticate, login
from django.contrib.auth.models import User
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from django.utils.decorators import method_decorator
from django.views.decorators.http import require_http_methods
from rest_framework import routers, status
from rest_framework.decorators import api_view
from rest_framework.response import Response
import json

from resolwe.api_urls import api_router as resolwe_router
from resolwe.flow.views import EntityViewSet

from resolwe_bio.filters import BioEntityFilter
from resolwe_bio.kb.views import FeatureViewSet, MappingSearchViewSet
from resolwe_bio.variants.views import (
    VariantAnnotationViewSet,
    VariantCallViewSet,
    VariantExperimentViewSet,
    VariantViewSet,
)

from .routers import SearchRouter


@api_view(['GET'])
def resdk_minimal_supported_version(request):
    """Return minimal supported ReSDK version."""
    return Response({"minimal_supported_version": "0.0.0"})


@api_view(['GET'])
def current_user(request):
    """Return current user info or anonymous user."""
    if hasattr(request, 'user') and request.user and hasattr(request.user, 'is_authenticated'):
        if request.user.is_authenticated:
            return Response({
                "id": getattr(request.user, 'id', 1),
                "username": getattr(request.user, 'username', 'anonymous'),
                "email": getattr(request.user, 'email', 'anonymous@localhost'),
                "first_name": getattr(request.user, 'first_name', ''),
                "last_name": getattr(request.user, 'last_name', ''),
            })
    
    # Return anonymous user
    return Response({
        "id": 1,
        "username": "anonymous", 
        "email": "anonymous@localhost",
        "first_name": "",
        "last_name": "",
    })


@api_view(['GET'])
def upload_config(request):
    """Return upload configuration."""
    return Response({
        "chunk_size": 1048576,  # 1MB chunks
        "upload_url": "/data/",
        "max_file_size": 1073741824,  # 1GB max file size
    })


@csrf_exempt
@require_http_methods(["POST"])
def upload_file(request):
    """Handle file upload endpoint."""
    try:
        # Get file from request
        if 'file' in request.FILES:
            uploaded_file = request.FILES['file']
            
            # Save file to uploads directory
            import os
            from django.conf import settings
            upload_dir = settings.FLOW_EXECUTOR["UPLOAD_DIR"]
            os.makedirs(upload_dir, exist_ok=True)
            
            file_path = os.path.join(upload_dir, uploaded_file.name)
            with open(file_path, 'wb+') as destination:
                for chunk in uploaded_file.chunks():
                    destination.write(chunk)
            
            return JsonResponse({
                "files": [{
                    "temp": uploaded_file.name,  # This is what resdk expects
                    "file_path": file_path,
                    "size": uploaded_file.size,
                }]
            })
        else:
            return JsonResponse({"error": "No file provided"}, status=400)
            
    except Exception as e:
        return JsonResponse({"error": str(e)}, status=500)


@csrf_exempt
@require_http_methods(["POST"])
def saml_auth_api_login(request):
    """Handle SAML authentication API login endpoint."""
    try:
        # Try JSON first, then form data
        if request.content_type == 'application/json':
            data = json.loads(request.body)
            username = data.get('username')
            password = data.get('password')
        else:
            username = request.POST.get('username') or request.POST.get('email')
            password = request.POST.get('password')
        
        # Create admin user if it doesn't exist
        if username == 'admin' and password == 'admin':
            user, created = User.objects.get_or_create(
                username='admin',
                defaults={
                    'email': 'admin@localhost',
                    'first_name': 'Admin',
                    'last_name': 'User',
                    'is_staff': True,
                    'is_superuser': True,
                }
            )
            if created:
                user.set_password('admin')
                user.save()
            
            # Authenticate and login
            user = authenticate(request, username='admin', password='admin')
            if user:
                login(request, user)
                
                # Create session if not exists
                if not request.session.session_key:
                    request.session.create()
                
                # Return response with cookies - resdk expects status 204 or 200
                response = JsonResponse({"status": "success"}, status=204)
                response.set_cookie('sessionid', request.session.session_key)
                # Generate a simple CSRF token
                from django.middleware.csrf import get_token
                csrf_token = get_token(request)
                response.set_cookie('csrftoken', csrf_token)
                return response
        
        return JsonResponse({"error": "Invalid credentials"}, status=401)
        
    except Exception as e:
        return JsonResponse({"error": str(e)}, status=400)

EntityViewSet.filterset_class = BioEntityFilter

api_router = routers.DefaultRouter(trailing_slash=False)
api_router.register(r"sample", EntityViewSet)
api_router.register(r"variant", VariantViewSet)
api_router.register(r"variant_annotations", VariantAnnotationViewSet)
api_router.register(r"variant_calls", VariantCallViewSet)
api_router.register(r"variant_experiment", VariantExperimentViewSet)

search_router = SearchRouter(trailing_slash=False)
search_router.register(r"kb/feature", FeatureViewSet, "kb_feature")
search_router.register(r"kb/mapping/search", MappingSearchViewSet, "kb_mapping_search")

urlpatterns = [
    path("api-auth/", include("rest_framework.urls", namespace="rest_framework")),
    # XXX: Temporary fix to work with Resolwe 2.0.0, which requires 'resolwe-api' namespace to be available when
    # reporting errors when running processes.
    path("api-resolwe/", include((resolwe_router.urls, "resolwe-api"))),
    # Add missing endpoints that resdk expects
    path("api/resdk_minimal_supported_version", resdk_minimal_supported_version, name="resdk_minimal_supported_version"),
    path("api/user", current_user, name="current_user"),
    path("api/upload_config", upload_config, name="upload_config"),
    # Add SAML authentication endpoint
    path("saml-auth/api-login/", saml_auth_api_login, name="saml_auth_api_login"),
    # Add upload endpoint
    path("upload/", upload_file, name="upload_file"),
    path(
        "api/",
        include(
            (
                api_router.urls + search_router.urls + resolwe_router.urls,
                "resolwebio-api",
            )
        ),
    ),
]
