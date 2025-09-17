#!/usr/bin/env bash
set -euo pipefail

# This script provisions a Debian/Ubuntu VPS with Docker (and compose plugin),
# Python 3, Node.js 20.x, and Caddy. It expects passwordless sudo for the user.

if ! command -v sudo >/dev/null 2>&1; then
  echo "sudo is required on the VPS" >&2
  exit 1
fi
if ! sudo -n true 2>/dev/null; then
  echo "Passwordless sudo required for provisioning. Configure NOPASSWD for user $(whoami)." >&2
  exit 1
fi

export DEBIAN_FRONTEND=noninteractive

if [ -f /etc/debian_version ] || [ -f /etc/lsb-release ]; then
  PKG_MGR=apt
else
  echo "Unsupported distribution (expecting Debian/Ubuntu)." >&2
  exit 1
fi

if [ "$PKG_MGR" = apt ]; then
  sudo apt-get update -y
  sudo apt-get install -y ca-certificates curl gnupg lsb-release

  # Docker Engine + Compose plugin
  if ! command -v docker >/dev/null 2>&1; then
    . /etc/os-release
    DIST_ID="$ID"
    DIST_CODENAME="$VERSION_CODENAME"
    sudo install -m 0755 -d /etc/apt/keyrings
    curl -fsSL https://download.docker.com/linux/${DIST_ID}/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
    sudo chmod a+r /etc/apt/keyrings/docker.gpg
    echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/${DIST_ID} ${DIST_CODENAME} stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
    sudo apt-get update -y
    sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
    sudo usermod -aG docker "$USER" || true
    sudo systemctl enable --now docker
  fi

  # Ensure Compose plugin is present even if Docker was preinstalled
  if ! docker compose version >/dev/null 2>&1; then
    . /etc/os-release
    DIST_ID="$ID"
    DIST_CODENAME="$VERSION_CODENAME"
    sudo install -m 0755 -d /etc/apt/keyrings || true
    curl -fsSL https://download.docker.com/linux/${DIST_ID}/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg || true
    sudo chmod a+r /etc/apt/keyrings/docker.gpg || true
    echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/${DIST_ID} ${DIST_CODENAME} stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null || true
    sudo apt-get update -y
    sudo apt-get install -y docker-buildx-plugin docker-compose-plugin
  fi

  # Python
  if ! command -v python3 >/dev/null 2>&1; then
    sudo apt-get install -y python3 python3-pip python3-venv
  else
    sudo apt-get install -y python3-pip python3-venv
  fi

  # Node.js (NodeSource 20.x LTS)
  if ! command -v node >/dev/null 2>&1; then
    curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash -
    sudo apt-get install -y nodejs
  fi

  # Caddy (reverse proxy)
  if ! command -v caddy >/dev/null 2>&1; then
    sudo apt-get install -y debian-keyring debian-archive-keyring apt-transport-https || true
    curl -fsSL https://dl.cloudsmith.io/public/caddy/stable/gpg.key | sudo gpg --dearmor -o /usr/share/keyrings/caddy-stable-archive-keyring.gpg
    echo "deb [signed-by=/usr/share/keyrings/caddy-stable-archive-keyring.gpg] https://dl.cloudsmith.io/public/caddy/stable/deb/debian any-version main" | sudo tee /etc/apt/sources.list.d/caddy-stable.list
    sudo apt-get update -y
    sudo apt-get install -y caddy
    sudo systemctl enable --now caddy || true
  fi
fi

echo "=== Versions ==="
docker --version || true
docker compose version || true
python3 --version || true
pip3 --version || true
node -v || true
npm -v || true
caddy version || true

echo "Provisioning complete. Reconnect may be required for docker group membership to take effect."


