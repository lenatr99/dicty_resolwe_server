#!/bin/bash
# Script to copy PostgreSQL Docker volume to VPS

set -euo pipefail

VPS_HOST="95.179.242.134"
VPS_USER="dictyapp"
LOCAL_VOLUME_NAME="backend_postgres_data"  # Local volume name
REMOTE_VOLUME_NAME="deploy_postgres_data"    # Remote volume name
BACKUP_FILE="/tmp/postgres_volume_backup.tar.gz"

# Check if local volume exists
if ! docker volume ls | grep -q "${LOCAL_VOLUME_NAME}"; then
    echo "Error: Local volume '${LOCAL_VOLUME_NAME}' not found!"
    echo "Available volumes:"
    docker volume ls
    exit 1
fi

echo "=== Creating volume backup ==="
# Create a compressed backup of the PostgreSQL volume
docker run --rm -v ${LOCAL_VOLUME_NAME}:/data -v /tmp:/backup alpine tar czf /backup/postgres_volume_backup.tar.gz -C /data .

echo "=== Checking backup size ==="
ls -lh ${BACKUP_FILE}

echo "=== Transferring volume backup to VPS ==="
# Transfer the volume backup to VPS
rsync -avz --progress ${BACKUP_FILE} ${VPS_USER}@${VPS_HOST}:/tmp/

echo "=== Restoring volume on VPS ==="
# Connect to VPS and restore the volume
ssh ${VPS_USER}@${VPS_HOST} << EOF
  # Stop services that use the database
  sudo systemctl stop dicty-backend || true
  cd ~/dicty_resolwe_server/deploy
  docker-compose -f docker-compose.db.yml down || true
  
  # Wait a moment for Docker to release resources
  sleep 5
  
  # Force remove any remaining containers using the volume
  docker ps -a --filter volume=${REMOTE_VOLUME_NAME} --format "{{.ID}}" | xargs -r docker rm -f || true
  
  # Prune system to clean up
  docker system prune -f || true
  
  # Try to remove the volume multiple times with delays
  for i in {1..5}; do
    if docker volume rm ${REMOTE_VOLUME_NAME} 2>/dev/null; then
      echo "Volume removed successfully"
      break
    else
      echo "Attempt \$i: Volume still in use, waiting..."
      sleep 3
    fi
  done
  
  # Create a new empty volume
  docker volume create ${REMOTE_VOLUME_NAME}
  
  # Restore the volume from backup
  docker run --rm -v ${REMOTE_VOLUME_NAME}:/data -v /tmp:/backup alpine tar xzf /backup/postgres_volume_backup.tar.gz -C /data
  
  # Start the database
  cd ~/dicty_resolwe_server/deploy
  docker-compose -f docker-compose.db.yml up -d
  
  # Wait for database to be ready
  echo "Waiting for database to start..."
  sleep 15
  
  # Check if database is responding
  for i in {1..30}; do
    if docker exec deploy-postgres-1 pg_isready -U resolwe >/dev/null 2>&1; then
      echo "Database is ready!"
      break
    else
      echo "Waiting for database... (\$i/30)"
      sleep 2
    fi
  done
  
  # Start the backend service
  sudo systemctl start dicty-backend
  
  echo "Volume restore complete!"
EOF

echo "=== Cleanup ==="
# Remove local backup file
rm -f ${BACKUP_FILE}

echo "Done! PostgreSQL volume has been transferred and restored on VPS."