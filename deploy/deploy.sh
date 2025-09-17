#!/bin/bash

# Dicty Resolwe Server Deployment Script
# This script deploys the application to the VPS server

set -e  # Exit on any error

# Configuration
PROJECT_DIR="/home/dictyapp/dicty_resolwe_server"
BACKEND_DIR="$PROJECT_DIR/backend"
FRONTEND_DIR="$PROJECT_DIR/frontend"
SERVICE_NAME="dicty-app"

echo "=== Starting Dicty Resolwe Server Deployment ==="
echo "Date: $(date)"
echo "User: $(whoami)"

# Update source code
echo "=== Updating source code ==="
cd $PROJECT_DIR

# Check if git repo exists, if not clone it
if [ ! -d ".git" ]; then
    echo "Git repository not found. Cloning..."
    cd /home/dictyapp
    git clone https://github.com/ltrnavove/dicty_resolwe_server.git
    cd dicty_resolwe_server
else
    echo "Git repository found. Pulling latest changes..."
    git fetch origin
    git reset --hard origin/master
fi

echo "Current commit: $(git rev-parse HEAD)"

# Stop existing services
echo "=== Stopping existing services ==="
sudo systemctl stop $SERVICE_NAME || true
sudo systemctl stop nginx || true

# Setup database
echo "=== Setting up database ==="
cd $BACKEND_DIR

# Setup environment file
if [ ! -f "$BACKEND_DIR/.env" ]; then
    echo "Setting up production environment file..."
    cp $PROJECT_DIR/deploy/production-env-template.txt $BACKEND_DIR/.env
    
    # Generate a random secret key
    SECRET_KEY=$(python3 -c 'from django.core.management.utils import get_random_secret_key; print(get_random_secret_key())')
    sed -i "s/your-super-secret-key-change-this-in-production-make-it-very-long-and-random/$SECRET_KEY/g" $BACKEND_DIR/.env
    
    echo "✅ Environment file created. Please review and update $BACKEND_DIR/.env with your specific settings."
fi

# Start Docker services for database using production compose
docker compose -f $PROJECT_DIR/deploy/docker-compose.prod.yml down || true
docker compose -f $PROJECT_DIR/deploy/docker-compose.prod.yml up -d

echo "Waiting for database to be ready..."
# Wait for PostgreSQL to be ready
until docker exec dicty-postgres pg_isready -U resolwe > /dev/null 2>&1; do
    echo "Waiting for database..."
    sleep 2
done

echo "Database is ready!"

# Setup Python environment for backend
echo "=== Setting up backend environment ==="
cd $BACKEND_DIR

# Create virtual environment if it doesn't exist
if [ ! -d "venv" ]; then
    python3 -m venv venv
fi

# Activate virtual environment and install dependencies
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt

# Load environment variables
if [ -f ".env" ]; then
    export $(cat .env | grep -v '^#' | xargs)
fi

# Database migrations and setup
echo "=== Running database migrations ==="
cd $BACKEND_DIR/resolwe_server

# Run Django migrations
python manage.py migrate

# Create superuser (skip if exists)
echo "=== Creating superuser ==="
python manage.py shell << 'EOF'
from django.contrib.auth import get_user_model
User = get_user_model()
if not User.objects.filter(username='admin').exists():
    User.objects.create_superuser('admin', 'admin@admin.com', 'admin123')
    print("Superuser created")
else:
    print("Superuser already exists")
EOF

# Register and prepare runtime
echo "=== Preparing runtime ==="
python manage.py register || true
python manage.py prepare_runtime || true
python manage.py collecttools || true

# Insert features if they don't exist
echo "=== Setting up features ==="
python manage.py shell << 'EOF'
from django.core.management import call_command
try:
    call_command('insert_features', 'dicty_features.tab')
    print("Features inserted")
except Exception as e:
    print(f"Features already exist or error: {e}")
EOF

# Setup annotations
echo "=== Setting up annotations ==="
cd $BACKEND_DIR
python setup_annotations.py || true

# Migrate database if needed
echo "=== Running database migration script ==="
if [ -f "reload_database.sh" ]; then
    # Update the script to use new VPS IP
    sed -i 's/91\.98\.119\.63/95.179.242.134/g' reload_database.sh
    sed -i 's/root/dictyapp/g' reload_database.sh
    echo "Database migration script updated for new VPS"
fi

# Build frontend
echo "=== Building frontend ==="
cd $FRONTEND_DIR

# Install Node.js and Yarn if not present
if ! command -v node &> /dev/null; then
    echo "Installing Node.js..."
    curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash -
    sudo apt-get install -y nodejs
fi

if ! command -v yarn &> /dev/null; then
    echo "Installing Yarn..."
    sudo npm install -g yarn
fi

# Install dependencies and build
echo "Installing frontend dependencies..."
yarn install --frozen-lockfile

echo "Building frontend..."
yarn build

# Setup Nginx configuration
echo "=== Setting up Nginx ==="
sudo tee /etc/nginx/sites-available/dicty-app << 'EOF'
server {
    listen 80;
    server_name 95.179.242.134;

    # Security headers
    add_header X-Frame-Options "SAMEORIGIN" always;
    add_header X-XSS-Protection "1; mode=block" always;
    add_header X-Content-Type-Options "nosniff" always;
    add_header Referrer-Policy "no-referrer-when-downgrade" always;
    add_header Content-Security-Policy "default-src 'self' http: https: data: blob: 'unsafe-inline'" always;

    # Serve frontend static files
    location / {
        root /home/dictyapp/dicty_resolwe_server/frontend/build;
        try_files $uri $uri/ /index.html;
        
        # Cache static assets
        location ~* \.(js|css|png|jpg|jpeg|gif|ico|svg)$ {
            expires 1y;
            add_header Cache-Control "public, immutable";
        }
    }

    # Proxy API requests to backend
    location /api/ {
        proxy_pass http://127.0.0.1:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        
        # Increase timeout for long-running requests
        proxy_read_timeout 300s;
        proxy_connect_timeout 75s;
    }

    # Serve Django admin
    location /admin/ {
        proxy_pass http://127.0.0.1:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    # Serve Django static files
    location /static/ {
        alias /home/dictyapp/dicty_resolwe_server/backend/staticfiles/;
        expires 1y;
        add_header Cache-Control "public";
    }

    # WebSocket support for Django Channels (if needed)
    location /ws/ {
        proxy_pass http://127.0.0.1:8000;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    # Security - deny access to sensitive files
    location ~ /\. {
        deny all;
    }
    
    location ~ /(\.git|node_modules|venv) {
        deny all;
    }
}
EOF

# Enable the site
sudo ln -sf /etc/nginx/sites-available/dicty-app /etc/nginx/sites-enabled/
sudo rm -f /etc/nginx/sites-enabled/default

# Test nginx configuration
sudo nginx -t

# Create systemd service for the Django app
echo "=== Creating systemd service ==="
sudo tee /etc/systemd/system/$SERVICE_NAME.service << EOF
[Unit]
Description=Dicty Resolwe Server Django Application
After=network.target
Requires=docker.service

[Service]
Type=simple
User=dictyapp
Group=dictyapp
WorkingDirectory=$BACKEND_DIR/resolwe_server
Environment=PATH=$BACKEND_DIR/venv/bin
ExecStart=$BACKEND_DIR/venv/bin/python manage.py runserver 0.0.0.0:8000
Restart=always
RestartSec=3
StandardOutput=journal
StandardError=journal

[Install]
WantedBy=multi-user.target
EOF

# Reload systemd and start services
echo "=== Starting services ==="
sudo systemctl daemon-reload
sudo systemctl enable $SERVICE_NAME
sudo systemctl start $SERVICE_NAME
sudo systemctl enable nginx
sudo systemctl start nginx

# Collect static files
echo "=== Collecting static files ==="
cd $BACKEND_DIR/resolwe_server
source ../venv/bin/activate
python manage.py collectstatic --noinput

# Wait for services to start
echo "=== Waiting for services to start ==="
sleep 10

# Check service status
echo "=== Checking service status ==="
sudo systemctl status $SERVICE_NAME --no-pager -l
sudo systemctl status nginx --no-pager -l

# Check if application is responding
echo "=== Testing application ==="
if curl -f http://localhost:8000/api/ > /dev/null 2>&1; then
    echo "✅ Backend API is responding"
else
    echo "❌ Backend API is not responding"
fi

if curl -f http://localhost/ > /dev/null 2>&1; then
    echo "✅ Frontend is responding"
else
    echo "❌ Frontend is not responding"
fi

echo "=== Deployment completed ==="
echo "Application should be available at: http://95.179.242.134"
echo "Django admin available at: http://95.179.242.134/admin"
echo ""
echo "Check logs with:"
echo "  sudo journalctl -u $SERVICE_NAME -f"
echo "  sudo tail -f /var/log/nginx/error.log"
