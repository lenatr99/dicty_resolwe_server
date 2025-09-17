#!/usr/bin/env bash
set -euo pipefail

# Backend deployment script - installs deps, runs migrations, starts service
# Usage: backend_deploy.sh [backend_dir] [env_file]

BACKEND_DIR=${1:-/home/$USER/dicty_resolwe_server/backend}
ENV_FILE=${2:-/home/$USER/dicty_env/backend.env}
VENV_DIR="$BACKEND_DIR/venv"

echo "=== Backend deployment starting ==="
echo "Backend dir: $BACKEND_DIR"
echo "Env file: $ENV_FILE"
echo "Venv dir: $VENV_DIR"

# Create backend env file directory if needed
mkdir -p "$(dirname "$ENV_FILE")"
touch "$ENV_FILE"

# Ensure we're in the backend directory
cd "$BACKEND_DIR"

# Create/update virtual environment
if [ ! -d "$VENV_DIR" ]; then
  echo "Creating Python virtual environment..."
  python3 -m venv "$VENV_DIR"
fi

# Activate virtual environment
source "$VENV_DIR/bin/activate"

# Upgrade pip and install dependencies
echo "Installing Python dependencies..."
pip install --upgrade pip
pip install -r requirements.txt
pip install gunicorn  # Production WSGI server

# Load environment variables
set -a
source "$ENV_FILE" 2>/dev/null || true
set +a

# Django setup
cd resolwe_server

# Collect static files
echo "Collecting static files..."
python manage.py collectstatic --noinput --clear

# Run database migrations
echo "Running database migrations..."
python manage.py migrate --noinput

# Create admin user if specified and doesn't exist
if [ -n "${ADMIN_USERNAME:-}" ] && [ -n "${ADMIN_PASSWORD:-}" ]; then
  echo "Ensuring admin user exists..."
  python manage.py shell -c "
from django.contrib.auth.models import User
username = '${ADMIN_USERNAME}'
password = '${ADMIN_PASSWORD}'
if not User.objects.filter(username=username).exists():
    User.objects.create_superuser(username, '${ADMIN_USERNAME}@localhost', password)
    print(f'Created admin user: {username}')
else:
    print(f'Admin user already exists: {username}')
" || true
fi

# Restart the backend service
echo "Restarting backend service..."
sudo systemctl daemon-reload
sudo systemctl restart dicty-backend || echo "Service not yet installed"
sudo systemctl enable dicty-backend || echo "Service not yet installed"

echo "=== Backend deployment complete ==="
