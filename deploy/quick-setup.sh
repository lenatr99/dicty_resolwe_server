#!/bin/bash

# Quick setup script for VPS - run this first to get the repository and setup
# This is a simplified version that you can copy-paste directly into the VPS

echo "=== Quick VPS Setup for Dicty Resolwe Server ==="
echo "VPS: 95.179.242.134"
echo "User: dictyapp"
echo ""

# Update system
echo "=== Updating system ==="
sudo apt update && sudo apt upgrade -y

# Install essential tools
echo "=== Installing essential tools ==="
sudo apt install -y git curl wget python3 python3-pip python3-venv docker.io docker-compose nginx ufw

# Add user to docker group
sudo usermod -aG docker $USER

# Clone repository
echo "=== Cloning repository ==="
cd /home/dictyapp
if [ -d "dicty_resolwe_server" ]; then
    echo "Repository already exists, updating..."
    cd dicty_resolwe_server
    git pull origin main
else
    git clone https://github.com/lenatr99/dicty_resolwe_server.git
    cd dicty_resolwe_server
fi

# Make scripts executable
chmod +x deploy/*.sh

echo "=== Repository cloned successfully ==="
echo ""
echo "Next steps:"
echo "1. Log out and log back in for docker group to take effect"
echo "2. Run: cd dicty_resolwe_server && ./deploy/server-setup.sh"
echo "3. Then run: ./deploy/deploy.sh"
echo ""
echo "Repository location: /home/dictyapp/dicty_resolwe_server"
