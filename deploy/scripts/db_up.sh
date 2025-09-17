#!/usr/bin/env bash
set -euo pipefail

COMPOSE_FILE=${1:-/home/$USER/dicty_resolwe_server/deploy/docker-compose.db.yml}
ENV_FILE=${2:-/home/$USER/dicty_resolwe_server/deploy/db.env}

# Ensure docker is available
if ! command -v docker >/dev/null 2>&1; then
  echo "Docker is required." >&2
  exit 1
fi

# Free port 5432 if another Docker container is using it
conflict_id=$(docker ps --format '{{.ID}} {{.Names}} {{.Ports}}' | awk '/127\.0\.0\.1:5432->|0\.0\.0\.0:5432->/ {print $1; exit}') || true
if [ -n "${conflict_id:-}" ]; then
  echo "Detected a container using host port 5432. Stopping container $conflict_id to free the port..."
  docker stop "$conflict_id" || true
  # Give Docker a moment to release the port
  sleep 2
fi

# Ensure compose plugin works; fallback to standalone compose v2 if needed
if ! docker compose version >/dev/null 2>&1; then
  echo "docker compose plugin not available; installing standalone docker-compose v2..."
  sudo curl -L "https://github.com/docker/compose/releases/download/v2.29.7/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
  sudo chmod +x /usr/local/bin/docker-compose
  comp_cmd="docker-compose"
else
  comp_cmd="docker compose"
fi

mkdir -p "$(dirname "$ENV_FILE")"
touch "$ENV_FILE"

echo "Bringing up PostgreSQL using $COMPOSE_FILE"
set -a
source "$ENV_FILE" || true
set +a

$comp_cmd --env-file "$ENV_FILE" -f "$COMPOSE_FILE" up -d

echo "Waiting for database to be ready..."
for i in {1..30}; do
  cname=$(docker ps --format '{{.Names}}' | grep '^.*postgres.*$' | head -n1 || true)
  if [ -n "$cname" ] && docker exec "$cname" pg_isready -U "${DB_USER:-resolwe}" >/dev/null 2>&1; then
    echo "Database is ready!"
    exit 0
  fi
  echo "Waiting ($i/30)..."
  sleep 2
done

echo "Database did not become ready in time." >&2
exit 1


