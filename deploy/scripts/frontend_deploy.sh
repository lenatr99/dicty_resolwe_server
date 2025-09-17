#!/usr/bin/env bash
set -euo pipefail

# Frontend deployment script - builds React app and deploys to nginx/caddy
# Usage: frontend_deploy.sh [frontend_dir] [build_output_dir]

FRONTEND_DIR=${1:-/home/$USER/dicty_resolwe_server/frontend}
BUILD_OUTPUT_DIR=${2:-/home/$USER/www/dicty}

echo "=== Frontend deployment starting ==="
echo "Frontend dir: $FRONTEND_DIR"
echo "Build output: $BUILD_OUTPUT_DIR"

# Ensure we're in the frontend directory
cd "$FRONTEND_DIR"

# Check if Node.js and npm are available
if ! command -v node >/dev/null 2>&1; then
  echo "Node.js is required but not found" >&2
  exit 1
fi

if ! command -v npm >/dev/null 2>&1; then
  echo "npm is required but not found" >&2
  exit 1
fi

echo "Node version: $(node -v)"
echo "npm version: $(npm -v)"

# Install dependencies
echo "Installing frontend dependencies..."
npm ci --only=production

# Build the application
echo "Building frontend for production..."
npm run build

# Create output directory
sudo mkdir -p "$BUILD_OUTPUT_DIR"

# Copy built files to web directory
echo "Deploying built files to $BUILD_OUTPUT_DIR..."
sudo cp -r build/* "$BUILD_OUTPUT_DIR"/
sudo chown -R dictyapp:dictyapp "$BUILD_OUTPUT_DIR"
sudo chmod -R 755 "$BUILD_OUTPUT_DIR"

# Create a simple index check
if [ -f "$BUILD_OUTPUT_DIR/index.html" ]; then
  echo "✅ Frontend deployment successful - index.html found"
  echo "Build size: $(du -sh "$BUILD_OUTPUT_DIR")"
else
  echo "❌ Frontend deployment failed - index.html not found"
  exit 1
fi

echo "=== Frontend deployment complete ==="
