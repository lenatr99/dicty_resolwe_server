# Dicty Resolwe Server

## Setup Instructions

### 1. Start Docker Services
```bash
cd backend && docker compose up
```

### 2. Initialize Database (new terminal)
```bash
cd backend/resolwe_server
python manage.py migrate
python manage.py createsuperuser --username admin --email admin@admin.com
python manage.py register
python manage.py prepare_runtime
python manage.py collecttools
python manage.py insert_features dicty_features.tab
python manage.py runserver
```

### 3. Start Workers (new terminal)
```bash
# Clean up existing processes
pkill -f runworker || true
pkill -f runlistener || true

# Start worker
cd backend/resolwe_server
python manage.py runworker resolwe-server.manager.control
```

### 4. Start Listener (new terminal)
```bash
cd backend/resolwe_server
python manage.py runlistener
```

### 5. Setup Annotations (new terminal)
```bash
cd backend
python setup_annotations.py
```