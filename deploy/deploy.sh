#!/usr/bin/env bash
set -euo pipefail

# Idempotent deploy script for VPS
# - Updates repo
# - Ensures Python 3.11 + venv, installs backend deps, runs migrations
# - Builds frontend and syncs to /var/www/dicty
# - Ensures Caddy is installed and configured
# - Ensures and restarts systemd service for backend (daphne)

PROJECT_DIR="${PROJECT_DIR:-$HOME/dicty_resolwe_server}"
BACKEND_DIR="$PROJECT_DIR/backend"
DJANGO_DIR="$BACKEND_DIR/resolwe_server"
VENV_DIR="$BACKEND_DIR/.venv"
FRONTEND_DIR="$PROJECT_DIR/frontend"
WEB_ROOT="/var/www/dicty"
SERVICE_NAME="dicty-backend"

echo "[deploy] Using PROJECT_DIR=$PROJECT_DIR"

update_repo() {
  # Repo content is delivered by rsync from GitHub Actions.
  # To avoid git ownership/safe.directory issues on the VPS, skip git here.
  # Ensure the project directory exists and proceed with build/deploy.
  mkdir -p "$PROJECT_DIR"
  echo "[deploy] Using rsync-delivered sources; skipping git operations"
}

ensure_python() {
  if ! command -v python3.11 >/dev/null 2>&1; then
    echo "[deploy] Installing Python 3.11"
    sudo apt update
    sudo apt install -y software-properties-common curl ca-certificates gnupg
    sudo add-apt-repository -y ppa:deadsnakes/ppa
    sudo apt update
    sudo apt install -y python3.11 python3.11-venv
  fi
}

ensure_swap() {
  # Create 4G swap if none exists (helps Node build)
  if ! grep -q "swapfile" /etc/fstab 2>/dev/null && [ ! -f /swapfile ]; then
    echo "[deploy] Creating 4G swapfile"
    sudo fallocate -l 4G /swapfile || sudo dd if=/dev/zero of=/swapfile bs=1M count=4096
    sudo chmod 600 /swapfile
    sudo mkswap /swapfile
    sudo swapon /swapfile
    echo '/swapfile none swap sw 0 0' | sudo tee -a /etc/fstab >/dev/null
  fi
}

setup_backend() {
  echo "[deploy] Backend setup"
  python3.11 -m venv "$VENV_DIR" || true
  source "$VENV_DIR/bin/activate"
  pip install -U pip setuptools wheel
  pip install -r "$BACKEND_DIR/requirements.txt"
  # Ensure .env exists at backend/.env (preferred path)
  if [ -f "$DJANGO_DIR/.env" ] && [ ! -f "$BACKEND_DIR/.env" ]; then
    mv "$DJANGO_DIR/.env" "$BACKEND_DIR/.env"
  fi
  # Migrate (safe)
  (cd "$DJANGO_DIR" && python manage.py migrate)
}

ensure_listener_keys() {
  # Generate ZMQ curve keys if not present in backend/.env
  local env_file="$BACKEND_DIR/.env"
  if [ ! -f "$env_file" ]; then
    touch "$env_file"
  fi
  if ! grep -q '^LISTENER_PUBLIC_KEY=' "$env_file" || ! grep -q '^LISTENER_PRIVATE_KEY=' "$env_file"; then
    echo "[deploy] Generating listener keys"
    local out
    out=$($VENV_DIR/bin/python - <<'PY'
import zmq
pub, priv = zmq.curve_keypair()
print(f"LISTENER_PUBLIC_KEY={pub.decode()}")
print(f"LISTENER_PRIVATE_KEY={priv.decode()}")
PY
)
    printf "%s\n" "$out" | sudo tee -a "$env_file" >/dev/null
  fi
}

ensure_node() {
  if ! command -v node >/dev/null 2>&1 || ! command -v npm >/dev/null 2>&1; then
    echo "[deploy] Installing Node.js 20"
    curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash -
    sudo apt install -y nodejs
  fi
  # Enable corepack/yarn
  if ! command -v yarn >/dev/null 2>&1; then
    corepack enable || true
  fi
}

build_frontend() {
  echo "[deploy] Building frontend"
  ensure_swap
  ensure_node
  export NODE_OPTIONS=${NODE_OPTIONS:---max-old-space-size=6144}
  (cd "$FRONTEND_DIR" && yarn install && yarn build)
  sudo mkdir -p "$WEB_ROOT"
  sudo rsync -ah --delete "$FRONTEND_DIR/build/" "$WEB_ROOT/"
}

ensure_caddy() {
  if ! command -v caddy >/dev/null 2>&1; then
    echo "[deploy] Installing Caddy"
    sudo apt update
    sudo apt install -y debian-keyring debian-archive-keyring apt-transport-https curl
    curl -1sLf 'https://dl.cloudsmith.io/public/caddy/stable/gpg.key' | sudo gpg --dearmor -o /usr/share/keyrings/caddy-stable-archive-keyring.gpg
    curl -1sLf 'https://dl.cloudsmith.io/public/caddy/stable/debian.deb.txt' | sudo tee /etc/apt/sources.list.d/caddy-stable.list >/dev/null
    sudo apt update
    sudo apt install -y caddy
  fi
  if [ -f "$PROJECT_DIR/deploy/Caddyfile" ]; then
    sudo install -m 0644 "$PROJECT_DIR/deploy/Caddyfile" /etc/caddy/Caddyfile
    sudo systemctl enable --now caddy
    sudo systemctl restart caddy
  fi
}

ensure_backend_service() {
  echo "[deploy] Ensuring systemd service $SERVICE_NAME"
  # Create a unit if not present or if template changed
  if [ -f "$PROJECT_DIR/deploy/backend.service" ]; then
    sudo install -m 0644 "$PROJECT_DIR/deploy/backend.service" \
      "/etc/systemd/system/${SERVICE_NAME}.service"
  else
    # Fallback minimal unit
    sudo tee "/etc/systemd/system/${SERVICE_NAME}.service" >/dev/null <<UNIT
[Unit]
Description=Dicty Resolwe Backend (Daphne)
After=network.target docker.service

[Service]
Type=simple
WorkingDirectory=$DJANGO_DIR
EnvironmentFile=$DJANGO_DIR/.env
ExecStart=$VENV_DIR/bin/daphne -b 0.0.0.0 -p 8000 resolwe_server.asgi:application
Restart=always

[Install]
WantedBy=multi-user.target
UNIT
  fi
  sudo systemctl daemon-reload
  sudo systemctl enable "$SERVICE_NAME"
  sudo systemctl restart "$SERVICE_NAME"
}

main() {
  update_repo
  ensure_python
  setup_backend
  build_frontend
  ensure_caddy
  ensure_listener_keys
  ensure_backend_service
  echo "[deploy] Done. Backend service: $SERVICE_NAME; Web root: $WEB_ROOT"
}

main "$@"


