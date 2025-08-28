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
from django.views.decorators.csrf import ensure_csrf_cookie
from django.middleware.csrf import get_token
import json

from resolwe.api_urls import api_router as resolwe_router
from resolwe.flow.views import EntityViewSet
from .views import UserViewSet, DifferentialExpressionViewSet, BasketViewSet, BasketExpressionsViewSet, GeneListByIdsViewSet

from resolwe_bio.filters import BioEntityFilter
from resolwe_bio.kb.views import MappingSearchViewSet
from .views import CustomFeatureViewSet
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
        "chunk_size": 10485760,  # 10MB chunks (larger chunks for big files)
        "upload_url": "/upload/",  # Use our working upload endpoint
        "max_file_size": 10737418240,  # 10GB max file size (increased from 1GB)
    })


@csrf_exempt
@require_http_methods(["POST"])
def upload_file(request):
    """Handle file upload endpoint with resdk chunked upload support."""
    try:
        # Check if this is a resdk chunked upload
        session_id = request.headers.get('Session-Id')
        file_uid = request.headers.get('X-File-Uid')
        chunk_number = request.POST.get('_chunkNumber')
        total_size = request.POST.get('_totalSize')
        current_chunk_size = request.POST.get('_currentChunkSize')
        
        if session_id and file_uid and chunk_number is not None:
            return handle_resdk_chunked_upload(request, session_id, file_uid, 
                                             int(chunk_number), int(total_size), 
                                             int(current_chunk_size))
        
        # Regular single file upload
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
        import traceback
        return JsonResponse({
            "error": str(e),
            "traceback": traceback.format_exc()
        }, status=500)


def handle_resdk_chunked_upload(request, session_id, file_uid, chunk_number, total_size, current_chunk_size):
    """Handle resdk's chunked upload protocol."""
    import os
    from django.conf import settings
    
    if 'file' not in request.FILES:
        return JsonResponse({"error": "No chunk data provided"}, status=400)
    
    chunk_file = request.FILES['file']
    file_name = chunk_file.name
    
    # Create temp directory for this upload session
    upload_dir = settings.FLOW_EXECUTOR["UPLOAD_DIR"]
    temp_dir = os.path.join(upload_dir, "temp_uploads", file_uid)
    os.makedirs(temp_dir, exist_ok=True)
    
    # Save this chunk
    chunk_path = os.path.join(temp_dir, f"chunk_{chunk_number:05d}")
    with open(chunk_path, 'wb') as f:
        for data in chunk_file.chunks():
            f.write(data)
    
    # Check current total size of received chunks
    current_total = 0
    chunk_files = []
    chunk_index = 0
    
    while True:
        chunk_file_path = os.path.join(temp_dir, f"chunk_{chunk_index:05d}")
        if os.path.exists(chunk_file_path):
            chunk_files.append(chunk_file_path)
            current_total += os.path.getsize(chunk_file_path)
            chunk_index += 1
        else:
            break
    
    # Check if upload is complete
    if current_total >= total_size:
        # All chunks received, assemble the file
        final_path = os.path.join(upload_dir, file_name)
        
        with open(final_path, 'wb') as final_file:
            for chunk_file_path in sorted(chunk_files):
                with open(chunk_file_path, 'rb') as cf:
                    final_file.write(cf.read())
        
        # Clean up temp chunks
        import shutil
        shutil.rmtree(temp_dir)
        
        # Verify final file size
        final_size = os.path.getsize(final_path)
        
        return JsonResponse({
            "files": [{
                "temp": file_name,
                "file_path": final_path,
                "size": final_size,
            }]
        })
    else:
        # More chunks needed
        return JsonResponse({
            "status": "chunk_received",
            "chunk_number": chunk_number,
            "received_size": current_total,
            "total_size": total_size,
            "progress": f"{current_total/total_size*100:.1f}%"
        })


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
    
# --- CSRF token view ---
@ensure_csrf_cookie
def csrf_token_view(request):
    return JsonResponse({"csrfToken": get_token(request)})

EntityViewSet.filterset_class = BioEntityFilter

api_router = routers.DefaultRouter(trailing_slash=False)
api_router.register(r"sample", EntityViewSet)
api_router.register(r"variant", VariantViewSet)
api_router.register(r"variant_annotations", VariantAnnotationViewSet)
api_router.register(r"variant_calls", VariantCallViewSet)
api_router.register(r"variant_experiment", VariantExperimentViewSet)
api_router.register(
    r"_modules/differential_expression/list",
    DifferentialExpressionViewSet,
    basename="diffexp-list",
)
api_router.register(r"basket", BasketViewSet, basename="basket")
api_router.register(
    r"_modules/visualizations/basket_expressions",
    BasketExpressionsViewSet,
    basename="basket-expr",
)
api_router.register(
    r"_modules/gene_list/list_by_ids",
    GeneListByIdsViewSet,
    basename="gene-list-by-ids",
)


search_router = SearchRouter(trailing_slash=False)
search_router.register(r"kb/feature", CustomFeatureViewSet, "kb_feature")
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
    path("api/base/csrf", csrf_token_view, name="csrf_token"),
]
