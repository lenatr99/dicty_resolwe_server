from django.contrib.auth import get_user_model
from rest_framework import serializers
from resolwe.flow.serializers import ResolweBaseSerializer
from rest_framework.serializers import ValidationError

User = get_user_model()



class UserSerializer(ResolweBaseSerializer):
    """Serializer for User objects."""
    
    class Meta:
        """Serializer configuration."""
        
        model = User
        read_only_fields = ("id", "date_joined", "last_login", "is_superuser")
        update_protected_fields = ("username", "email")
        fields = read_only_fields + update_protected_fields + (
            "first_name", 
            "last_name", 
            "is_active", 
            "is_staff"
        )
        
    def create(self, validated_data):
        """Create a new user instance."""
        # Ensure only staff can create users
        request = self.context.get('request')
        if not request.user.is_staff:
            raise ValidationError("Only staff users can create new users.")
        
        return super().create(validated_data)



class BasketCreateSerializer(serializers.Serializer):
    # what the client actually sends
    annotation_version = serializers.CharField()
    only_existing_organism = serializers.BooleanField()
    samples = serializers.ListField(child=serializers.IntegerField(min_value=1))

    # nice-to-have fields for browsable API (optional)
    name = serializers.CharField(required=False, allow_blank=True, default="")
    read_only = serializers.BooleanField(required=False, default=False)
    contributor = serializers.SerializerMethodField(read_only=True)

    def get_contributor(self, obj):
        u = self.context["request"].user
        return {"first_name": u.first_name, "last_name": u.last_name, "username": u.username}
