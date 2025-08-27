"""
ASGI config for resolwe_server project.

It exposes the ASGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/5.1/howto/deployment/asgi/
"""

import os
from django.core.asgi import get_asgi_application
from channels.routing import ProtocolTypeRouter, URLRouter
from channels.auth import AuthMiddlewareStack

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'settings')

django_asgi_app = get_asgi_application()

# Import Resolwe flow components
from channels.routing import ChannelNameRouter
from resolwe.flow.managers import state
from resolwe.flow.managers.consumer import ManagerConsumer, HealtCheckConsumer, CHANNEL_HEALTH_CHECK

# Create custom channel routing with proper ASGI instantiation
custom_channel_routing = ProtocolTypeRouter({
    "channel": ChannelNameRouter({
        state.MANAGER_CONTROL_CHANNEL: ManagerConsumer.as_asgi(),
        CHANNEL_HEALTH_CHECK: HealtCheckConsumer.as_asgi(),
    }),
})

application = ProtocolTypeRouter({
    "http": django_asgi_app,
    "websocket": AuthMiddlewareStack(
        URLRouter([
            # Add WebSocket URL patterns here if needed
        ])
    ),
    "channel": custom_channel_routing.application_mapping["channel"],
})