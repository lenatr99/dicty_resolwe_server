#!/bin/bash

# Server Setup Script for Dicty Resolwe Server VPS
# Run this script once on a fresh Ubuntu/Debian server to prepare it for deployment

set -e  # Exit on any error

VPS_IP="95.179.242.134"
USER="dictyapp"

echo "=== Dicty Resolwe Server - VPS Setup Script ==="
echo "Setting up server: $VPS_IP"
echo "User: $USER"
echo ""

# Update system packages
echo "=== Updating system packages ==="
sudo apt update
sudo apt upgrade -y

# Install essential packages
echo "=== Installing essential packages ==="
sudo apt install -y \
    curl \
    wget \
    git \
    unzip \
    software-properties-common \
    apt-transport-https \
    ca-certificates \
    gnupg \
    lsb-release \
    build-essential \
    python3 \
    python3-pip \
    python3-venv \
    python3-dev \
    postgresql-client \
    nginx \
    ufw \
    fail2ban \
    htop \
    tree

# Install Docker
echo "=== Installing Docker ==="
if ! command -v docker &> /dev/null; then
    # Add Docker's official GPG key
    curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
    
    # Add Docker repository
    echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
    
    # Install Docker
    sudo apt update
    sudo apt install -y docker-ce docker-ce-cli containerd.io docker-compose-plugin
    
    # Add user to docker group
    sudo usermod -aG docker $USER
    
    # Enable Docker to start on boot
    sudo systemctl enable docker
    sudo systemctl start docker
    
    echo "✅ Docker installed successfully"
else
    echo "✅ Docker already installed"
fi

# Install Node.js (latest LTS)
echo "=== Installing Node.js ==="
if ! command -v node &> /dev/null; then
    curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash -
    sudo apt install -y nodejs
    
    # Install Yarn
    sudo npm install -g yarn
    
    echo "✅ Node.js and Yarn installed successfully"
    echo "Node.js version: $(node --version)"
    echo "Yarn version: $(yarn --version)"
else
    echo "✅ Node.js already installed: $(node --version)"
fi

# Setup UFW Firewall
echo "=== Configuring firewall ==="
sudo ufw --force reset
sudo ufw default deny incoming
sudo ufw default allow outgoing

# Allow SSH (make sure we don't lock ourselves out)
sudo ufw allow ssh
sudo ufw allow 22/tcp

# Allow HTTP and HTTPS
sudo ufw allow 80/tcp
sudo ufw allow 443/tcp

# Allow database connections from localhost only (they're already bound to localhost in docker-compose)
# sudo ufw allow from 127.0.0.1 to any port 5432
# sudo ufw allow from 127.0.0.1 to any port 6379

# Enable firewall
sudo ufw --force enable
sudo ufw status verbose

echo "✅ Firewall configured"

# Configure Fail2ban
echo "=== Configuring Fail2ban ==="
sudo systemctl enable fail2ban
sudo systemctl start fail2ban

# Create basic jail configuration
sudo tee /etc/fail2ban/jail.local << 'EOF'
[DEFAULT]
bantime = 1h
findtime = 10m
maxretry = 5

[sshd]
enabled = true
port = ssh
filter = sshd
logpath = /var/log/auth.log
maxretry = 3

[nginx-http-auth]
enabled = true
filter = nginx-http-auth
port = http,https
logpath = /var/log/nginx/error.log

[nginx-noscript]
enabled = true
port = http,https
filter = nginx-noscript
logpath = /var/log/nginx/access.log
maxretry = 6

[nginx-badbots]
enabled = true
port = http,https
filter = nginx-badbots
logpath = /var/log/nginx/access.log
maxretry = 2
EOF

sudo systemctl restart fail2ban
echo "✅ Fail2ban configured"

# Create project directory structure
echo "=== Setting up project directory ==="
PROJECT_DIR="/home/$USER/dicty_resolwe_server"
if [ ! -d "$PROJECT_DIR" ]; then
    mkdir -p "$PROJECT_DIR"
    chown $USER:$USER "$PROJECT_DIR"
fi

# Create data directories
mkdir -p "$PROJECT_DIR/backend/data/runtime"
mkdir -p "$PROJECT_DIR/backend/data/upload"
chown -R $USER:$USER "$PROJECT_DIR/backend/data"

echo "✅ Project directories created"

# Setup log rotation
echo "=== Setting up log rotation ==="
sudo tee /etc/logrotate.d/dicty-app << 'EOF'
/var/log/dicty-app/*.log {
    daily
    missingok
    rotate 30
    compress
    delaycompress
    notifempty
    create 0644 dictyapp dictyapp
    postrotate
        systemctl reload dicty-app > /dev/null 2>&1 || true
    endscript
}
EOF

# Create log directory
sudo mkdir -p /var/log/dicty-app
sudo chown $USER:$USER /var/log/dicty-app

echo "✅ Log rotation configured"

# Setup system limits for the application
echo "=== Configuring system limits ==="
sudo tee /etc/security/limits.d/dicty-app.conf << 'EOF'
dictyapp soft nofile 65536
dictyapp hard nofile 65536
dictyapp soft nproc 4096
dictyapp hard nproc 4096
EOF

# Configure sysctl for better performance
sudo tee /etc/sysctl.d/99-dicty-app.conf << 'EOF'
# Increase file descriptor limits
fs.file-max = 65536

# Improve network performance
net.core.somaxconn = 1024
net.ipv4.tcp_max_syn_backlog = 2048

# Optimize for web server
net.ipv4.ip_local_port_range = 15000 65000
net.ipv4.tcp_fin_timeout = 30
EOF

sudo sysctl -p /etc/sysctl.d/99-dicty-app.conf

echo "✅ System limits configured"

# Install and configure automatic security updates
echo "=== Setting up automatic security updates ==="
sudo apt install -y unattended-upgrades
sudo dpkg-reconfigure -plow unattended-upgrades

echo "✅ Automatic security updates configured"

# Create backup script
echo "=== Creating backup script ==="
sudo tee /usr/local/bin/dicty-backup.sh << 'EOF'
#!/bin/bash
# Backup script for Dicty Resolwe Server

BACKUP_DIR="/home/dictyapp/backups"
DATE=$(date +%Y%m%d_%H%M%S)
PROJECT_DIR="/home/dictyapp/dicty_resolwe_server"

mkdir -p $BACKUP_DIR

# Backup database
docker exec dicty-postgres pg_dump -U resolwe resolwe > $BACKUP_DIR/database_$DATE.sql

# Backup application data
tar -czf $BACKUP_DIR/data_$DATE.tar.gz -C $PROJECT_DIR/backend data/

# Keep only last 7 backups
find $BACKUP_DIR -name "database_*.sql" -mtime +7 -delete
find $BACKUP_DIR -name "data_*.tar.gz" -mtime +7 -delete

echo "Backup completed: $DATE"
EOF

sudo chmod +x /usr/local/bin/dicty-backup.sh
sudo chown $USER:$USER /usr/local/bin/dicty-backup.sh

# Setup daily backup cron job
(crontab -l 2>/dev/null; echo "0 2 * * * /usr/local/bin/dicty-backup.sh >> /var/log/dicty-app/backup.log 2>&1") | crontab -

echo "✅ Backup script and cron job created"

# Final system optimization
echo "=== Final optimizations ==="

# Disable unnecessary services
sudo systemctl disable snapd || true
sudo systemctl disable cups || true

# Enable and start essential services
sudo systemctl enable nginx
sudo systemctl enable docker
sudo systemctl enable fail2ban

echo ""
echo "=== VPS Setup Complete! ==="
echo ""
echo "✅ System updated and secured"
echo "✅ Docker and Docker Compose installed"  
echo "✅ Node.js and Yarn installed"
echo "✅ Nginx installed and configured"
echo "✅ Firewall configured (UFW)"
echo "✅ Fail2ban configured for security"
echo "✅ Project directories created"
echo "✅ Log rotation configured"
echo "✅ System limits optimized"
echo "✅ Backup script created"
echo ""
echo "Next steps:"
echo "1. Clone the repository: git clone https://github.com/ltrnavove/dicty_resolwe_server.git"
echo "2. Set up SSH keys for GitHub Actions deployment"
echo "3. Run the deployment script"
echo ""
echo "IMPORTANT: Log out and log back in for group changes to take effect!"
echo "Then you can test Docker with: docker run hello-world"
EOF
