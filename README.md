# Dicty Resolwe Server

## Setup Instructions

### 1. Start Docker Services
```bash
cd backend
docker compose up
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

### 3. Kill existing processes
```bash
pkill -f runworker || true
pkill -f runlistener || true
```

### 4. Start Listener (new terminal)
```bash
cd backend/resolwe_server
python manage.py runlistener
```

### 5. Start workers (new terminal)
```bash
cd backend/resolwe_server
python manage.py runworker resolwe-server.manager.control
```

### 6. Setup Annotations (new terminal)
```bash
cd backend
python setup_annotations.py
```

### 7. Add data
```bash
cd backend/resolwe_server
python test3.py
```

## Start frontend

```bash
cd frontend
yarn start
```

## Copy local data to backend
```bash
cd backend
./reload_database.sh
```