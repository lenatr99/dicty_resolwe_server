#!/usr/bin/env bash
set -euo pipefail

COMPOSE_FILE=${1:-/home/$USER/dicty_resolwe_server/deploy/docker-compose.db.yml}
ENV_FILE=${2:-/home/$USER/dicty_resolwe_server/deploy/db.env}

# Ensure docker is available
if ! command -v docker >/dev/null 2>&1; then
  echo "Docker is required." >&2
  exit 1
fi

# Free port 5432 if in use by Docker or host service
# 1) Stop any docker containers publishing 5432
conflict_ids=$(docker ps --filter "publish=5432" -q || true)
if [ -n "${conflict_ids:-}" ]; then
  echo "Stopping containers publishing host port 5432: $conflict_ids"
  docker stop $conflict_ids || true
  sleep 2
fi

# 2) If still in use by host process (e.g., system PostgreSQL), stop it gracefully
if command -v ss >/dev/null 2>&1; then
  in_use=$(sudo ss -ltnp 2>/dev/null | awk '$4 ~ /:5432$/ {print $0; exit}') || true
else
  in_use=$(sudo lsof -iTCP:5432 -sTCP:LISTEN -Pn 2>/dev/null | sed -n '2p') || true
fi

if [ -n "${in_use:-}" ]; then
  echo "Port 5432 is in use by host process:"
  echo "$in_use"
  echo "Attempting to stop system PostgreSQL service if running..."
  sudo systemctl stop postgresql || sudo systemctl stop postgresql@16-main || sudo systemctl stop postgresql-16 || true
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


