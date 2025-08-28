"""
Django settings for the new Resolwe server project.
"""

import os
from pathlib import Path
from decouple import config

from dotenv import load_dotenv


# Build paths inside the project like this: BASE_DIR / 'subdir'.
BASE_DIR = Path(__file__).resolve().parent.parent
load_dotenv(BASE_DIR / ".env")  # adjust path if your .env is elsewhere


# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = config("SECRET_KEY", default="your-secret-key-change-in-production")

# SECURITY WARNING: don't run with debug turned on in production!
# Set to False to use tools volume instead of conda environment paths
DEBUG = config("DEBUG", default=False, cast=bool)

ALLOWED_HOSTS = ["*"]

# Application definition
INSTALLED_APPS = [
    "django.contrib.admin",
    "django.contrib.auth",
    "django.contrib.contenttypes",
    "django.contrib.sessions",
    "django.contrib.messages",
    "django.contrib.staticfiles",
    # Third-party apps
    "rest_framework",
    "django_filters",
    "corsheaders",
    "channels",
    "guardian",
    # Resolwe apps
    "resolwe",
    "resolwe.permissions",
    "resolwe.flow",
    "resolwe.storage",
    "resolwe.observers",
    # Local apps
    "resolwe_server.base",
    "resolwe_bio",
    "resolwe_bio.kb",
    "resolwe_bio.variants",
]

FLOW_PROCESSES_RUNTIMES = (
    "resolwe.process.runtime.Process",
    "resolwe_bio.process.runtime.ProcessBio",
)



MIDDLEWARE = [
    "corsheaders.middleware.CorsMiddleware",
    "django.middleware.security.SecurityMiddleware",
    "django.contrib.sessions.middleware.SessionMiddleware",
    "django.middleware.common.CommonMiddleware",
    "django.middleware.csrf.CsrfViewMiddleware",
    "django.contrib.auth.middleware.AuthenticationMiddleware",
    "django.contrib.messages.middleware.MessageMiddleware",
    "django.middleware.clickjacking.XFrameOptionsMiddleware",
]

ROOT_URLCONF = "resolwe_server.urls"

TEMPLATES = [
    {
        "BACKEND": "django.template.backends.django.DjangoTemplates",
        "DIRS": [],
        "APP_DIRS": True,
        "OPTIONS": {
            "context_processors": [
                "django.template.context_processors.debug",
                "django.template.context_processors.request",
                "django.contrib.auth.context_processors.auth",
                "django.contrib.messages.context_processors.messages",
            ],
        },
    },
]

WSGI_APPLICATION = "resolwe_server.wsgi.application"
ASGI_APPLICATION = "resolwe_server.asgi.application"

# Database
DATABASES = {
    "default": {
        "ENGINE": "django.db.backends.postgresql",
        "NAME": config("DB_NAME", default="resolwe"),
        "USER": config("DB_USER", default="resolwe"),
        "PASSWORD": config("DB_PASSWORD", default="resolwe"),
        "HOST": config("DB_HOST", default="localhost"),
        "PORT": config("DB_PORT", default="5432"),
    }
}

# Password validation
AUTH_PASSWORD_VALIDATORS = [
    {
        "NAME": "django.contrib.auth.password_validation.UserAttributeSimilarityValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.MinimumLengthValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.CommonPasswordValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.NumericPasswordValidator",
    },
]

# Internationalization
LANGUAGE_CODE = "en-us"
TIME_ZONE = "UTC"
USE_I18N = True
USE_TZ = True

# Static files (CSS, JavaScript, Images)
STATIC_URL = "/static/"
STATIC_ROOT = BASE_DIR / "staticfiles"

# Media files
MEDIA_URL = "/media/"
MEDIA_ROOT = BASE_DIR / "media"

# Default primary key field type
DEFAULT_AUTO_FIELD = "django.db.models.BigAutoField"

# File upload settings for large files (like 3.6GB h5ad files)
FILE_UPLOAD_MAX_MEMORY_SIZE = 10 * 1024 * 1024 * 1024  # 10GB
DATA_UPLOAD_MAX_MEMORY_SIZE = 10 * 1024 * 1024 * 1024  # 10GB
FILE_UPLOAD_PERMISSIONS = 0o644

# CORS settings
CORS_ALLOW_ALL_ORIGINS = True  # Only for development
CORS_ALLOW_CREDENTIALS = True

# Redis configuration
REDIS_CONNECTION = {
    "host": config("REDIS_HOST", default="localhost"),
    "port": config("REDIS_PORT", default=6379, cast=int),
    "db": config("REDIS_DB", default=1, cast=int),
}

# Celery configuration
CELERY_BROKER_URL = f"redis://{REDIS_CONNECTION['host']}:{REDIS_CONNECTION['port']}/{REDIS_CONNECTION['db']}"
CELERY_RESULT_BACKEND = CELERY_BROKER_URL

# Channels configuration
CHANNEL_LAYERS = {
    "default": {
        "BACKEND": "channels_redis.core.RedisChannelLayer",
        "CONFIG": {
            "hosts": [(REDIS_CONNECTION["host"], REDIS_CONNECTION["port"])],
        },
    },
}

# Add the channel routing for Django Channels
CHANNEL_ROUTING = "resolwe_server.routing.channel_routing"

# Resolwe configuration - WORKING STANDALONE EXECUTION
# PROPER DOCKER CONTAINERIZATION SETUP
FLOW_MANAGER_DISABLE_AUTO_CALLS = False

# Docker network configuration for local execution
FLOW_DOCKER_NETWORK = "bridge"

FLOW_EXECUTOR = {
    "NAME": "resolwe.flow.executors.docker",  # Use Docker executor
    "DATA_DIR": str(BASE_DIR / "data" / "data"),
    "UPLOAD_DIR": str(BASE_DIR / "data" / "upload"),
    "RUNTIME_DIR": str(BASE_DIR / "data" / "runtime"),
    "LISTENER_CONNECTION": {
        "port": 53892,
        "protocol": "tcp",
    },
    "NETWORK": FLOW_DOCKER_NETWORK,
}

# Keep data after processing to ensure output files are transferred back
FLOW_MANAGER_KEEP_DATA = True

# Note: Dispatcher source code patched directly in resolwe/flow/managers/dispatcher.py

# Docker images required for Resolwe execution
FLOW_DOCKER_COMMUNICATOR_IMAGE = "public.ecr.aws/s4q6j6e8/resolwe/com:latest"
FLOW_DOCKER_DEFAULT_PROCESSING_CONTAINER_IMAGE = "public.ecr.aws/s4q6j6e8/resolwe/base:ubuntu-20.04"

# Volume configuration for Docker containers - must match storage settings format
FLOW_VOLUMES = {
    "runtime": {
        "type": "host_path",
        "config": {
            "path": str(BASE_DIR / "data" / "runtime"),
            "name": "runtime",
            "read_only": False,
        }
    },
    "processing": {
        "type": "host_path",
        "config": {
            "path": str(BASE_DIR / "data" / "data"),
            "name": "processing",
        }
    },
    "input": {
        "type": "host_path",
        "config": {
            "path": str(BASE_DIR / "data" / "data"),
            "name": "input",
            "read_only": False,
        }
    },
    "tools": {
        "type": "host_path",
        "config": {
            "path": str(BASE_DIR / "data" / "tools"),
            "name": "tools",
            "read_only": True,
        }
    },
    "secrets": {
        "type": "temporary_directory",
        "config": {
            "name": "secrets",
            "selinux_label": "z",
        }
    },
    "sockets": {
        "type": "temporary_directory",
        "config": {
            "name": "sockets",
            "selinux_label": "z",
        }
    },
}

# Communication container listener connection for local execution
COMMUNICATION_CONTAINER_LISTENER_CONNECTION = {
    "local": "127.0.0.1"
}

# Force non-DEBUG behavior for tools to use collected tools volume instead of conda paths
# With DEBUG=False, get_tools_paths() will automatically use the tools volume from FLOW_VOLUMES


FLOW_MANAGER = {
    "REDIS_PREFIX": "resolwe-server.manager",
    "REDIS_CONNECTION": REDIS_CONNECTION,
    "DISPATCHER_MAPPING": {
        "Interactive": "resolwe.flow.managers.workload_connectors.local",
        "Batch": "resolwe.flow.managers.workload_connectors.local",
    },
}

FLOW_API = {
    "PERMISSIONS": "resolwe.permissions.permissions",
}

# Enable auto calls to allow manager to process data objects (duplicate removed - defined earlier)

# Storage configuration
STORAGE_CONNECTORS = {
    "data": {
        "connector": "resolwe.storage.connectors.localconnector.LocalFilesystemConnector",
        "config": {
            "priority": 100,
            "path": str(BASE_DIR / "data" / "data"),
            "selinux_label": "z",
            "public_url": "local_data",
        },
    },
    "upload": {
        "connector": "resolwe.storage.connectors.localconnector.LocalFilesystemConnector",
        "config": {
            "priority": 100,
            "path": str(BASE_DIR / "data" / "upload"),
            "selinux_label": "z",
        },
    },
}

FLOW_STORAGE = {
    "data": {"connectors": ["data"]},
    "upload": {"connectors": ["upload"]},
    # "inputs": {"connectors": ["local"]},
}

# Storage manager configuration - ensure directories are created
FLOW_MANAGER = {
    "REDIS_PREFIX": "resolwe-server.manager",
    "REDIS_CONNECTION": REDIS_CONNECTION,
    "DISPATCHER_MAPPING": {
        "Interactive": "resolwe.flow.managers.workload_connectors.local",
        "Batch": "resolwe.flow.managers.workload_connectors.local",
    },
    "ENSURE_DATA_DIR": True,  # Force directory creation
}

# Ensure storage directories are created before processing
STORAGE_ENSURE_DIRS = True

FLOW_EXPRESSION_ENGINES = [
    {
        "ENGINE": "resolwe.flow.expression_engines.jinja",
        "CUSTOM_FILTERS": [
            "resolwe_bio.expression_filters.sample",
            "resolwe_bio.expression_filters.relation",
        ],
    },
]

FLOW_EXECUTION_ENGINES = [
    "resolwe.flow.execution_engines.bash",
    "resolwe.flow.execution_engines.workflow",
    "resolwe.flow.execution_engines.python",
]

FLOW_PROCESSES_FINDERS = (
    "resolwe.flow.finders.FileSystemProcessesFinder",
    "resolwe.flow.finders.AppDirectoriesFinder",
)


# Django REST Framework
REST_FRAMEWORK = {
    "DEFAULT_AUTHENTICATION_CLASSES": (
        "rest_framework.authentication.SessionAuthentication",
    ),
    "DEFAULT_FILTER_BACKENDS": (
        "django_filters.rest_framework.DjangoFilterBackend",
        # 'resolwe.permissions.filters.ResolwePermissionsFilter',
    ),
    "DEFAULT_PAGINATION_CLASS": "rest_framework.pagination.LimitOffsetPagination",
    "PAGE_SIZE": 100,
    "TEST_REQUEST_DEFAULT_FORMAT": "json",
}

# Guardian
AUTHENTICATION_BACKENDS = (
    "django.contrib.auth.backends.ModelBackend",
    "guardian.backends.ObjectPermissionBackend",
)

ANONYMOUS_USER_NAME = "AnonymousUser"

# Logging
LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "standard": {
            "format": "%(asctime)s - %(levelname)s - %(name)s[%(process)s]: %(message)s",
        },
    },
    "handlers": {
        "console": {
            "class": "logging.StreamHandler",
            "level": "INFO",
            "formatter": "standard",
        },
    },
    "root": {
        "handlers": ["console"],
        "level": "INFO",
    },
    "loggers": {
        "resolwe": {
            "handlers": ["console"],
            "level": "DEBUG",
            "propagate": False,
        },
    },
}
