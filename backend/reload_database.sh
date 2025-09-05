#!/bin/bash
# Script to copy PostgreSQL Docker volume to VPS

VPS_HOST="91.98.119.63"
VPS_USER="root"
VOLUME_NAME="backend_postgres_data"
BACKUP_FILE="/tmp/postgres_volume_backup.tar.gz"

echo "=== Creating volume backup ==="
# Create a compressed backup of the PostgreSQL volume
docker run --rm -v ${VOLUME_NAME}:/data -v /tmp:/backup alpine tar czf /backup/postgres_volume_backup.tar.gz -C /data .

echo "=== Checking backup size ==="
ls -lh ${BACKUP_FILE}

echo "=== Transferring volume backup to VPS ==="
# Transfer the volume backup to VPS
rsync -avz --progress ${BACKUP_FILE} ${VPS_USER}@${VPS_HOST}:/tmp/

echo "=== Restoring volume on VPS ==="
# Connect to VPS and restore the volume
ssh ${VPS_USER}@${VPS_HOST} << 'EOF'
  # Stop services that use the database
  systemctl stop dicty-backend
  docker compose -f ~/dicty_resolwe_server/backend/docker-compose.yml down
  
  # Wait a moment for Docker to release resources
  sleep 5
  
  # Force remove any remaining containers using the volume
  docker ps -a --filter volume=backend_postgres_data --format "{{.ID}}" | xargs -r docker rm -f
  
  # Prune system to clean up
  docker system prune -f
  
  # Try to remove the volume multiple times with delays
  for i in {1..5}; do
    if docker volume rm backend_postgres_data 2>/dev/null; then
      echo "Volume removed successfully"
      break
    else
      echo "Attempt $i: Volume still in use, waiting..."
      sleep 3
    fi
  done
  
  # Create a new empty volume
  docker volume create backend_postgres_data
  
  # Restore the volume from backup
  docker run --rm -v backend_postgres_data:/data -v /tmp:/backup alpine tar xzf /backup/postgres_volume_backup.tar.gz -C /data
  
  # Start the database
  cd ~/dicty_resolwe_server/backend
  docker compose up -d
  
  # Wait for database to be ready
  echo "Waiting for database to start..."
  sleep 15
  
  # Check if database is responding
  for i in {1..30}; do
    if docker exec backend-postgres-1 pg_isready -U resolwe >/dev/null 2>&1; then
      echo "Database is ready!"
      break
    else
      echo "Waiting for database... ($i/30)"
      sleep 2
    fi
  done
  
  # Start the backend service
  systemctl start dicty-backend
  
  echo "Volume restore complete!"
EOF

echo "=== Cleanup ==="
# Remove local backup file
rm -f ${BACKUP_FILE}

echo "Done! PostgreSQL volume has been transferred and restored on VPS."