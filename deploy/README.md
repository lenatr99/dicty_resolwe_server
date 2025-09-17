# Dicty Resolwe Server Deployment Guide

This directory contains all the deployment scripts and configurations for the Dicty Resolwe Server on VPS `95.179.242.134`.

## Quick Start

### 1. Initial VPS Setup (Run once)

Connect to your VPS as root and create the dictyapp user:
```bash
ssh root@95.179.242.134
adduser dictyapp
usermod -aG sudo dictyapp
su - dictyapp
```

Run the server setup script to install all dependencies:
```bash
wget https://raw.githubusercontent.com/ltrnavove/dicty_resolwe_server/master/deploy/server-setup.sh
chmod +x server-setup.sh
./server-setup.sh
```

### 2. Setup SSH Keys for GitHub Actions

Generate SSH key pair on your local machine:
```bash
ssh-keygen -t ed25519 -C "github-actions-dicty" -f ~/.ssh/dicty_deploy_key
```

Copy the public key to the VPS:
```bash
ssh-copy-id -i ~/.ssh/dicty_deploy_key.pub dictyapp@95.179.242.134
```

Add the private key to GitHub repository secrets:
1. Go to your GitHub repository
2. Settings → Secrets and variables → Actions
3. Add new secret named `VPS_SSH_KEY`
4. Paste the content of `~/.ssh/dicty_deploy_key` (private key)

### 3. Initial Manual Deployment

Clone the repository on the VPS:
```bash
ssh dictyapp@95.179.242.134
git clone https://github.com/ltrnavove/dicty_resolwe_server.git
cd dicty_resolwe_server
chmod +x deploy/deploy.sh
./deploy/deploy.sh
```

### 4. Automatic Deployment

Every push to the `master` branch will now trigger automatic deployment via GitHub Actions.

## Files Description

### Core Deployment Files
- `deploy.sh` - Main deployment script that handles the entire deployment process
- `docker-compose.prod.yml` - Production Docker Compose configuration for database services
- `production-env-template.txt` - Environment variables template for production

### Setup Files
- `server-setup.sh` - Initial VPS setup script (run once)
- `dicty-app.service` - Systemd service configuration
- `Caddyfile` - Alternative web server configuration (with automatic HTTPS)

### Monitoring & Maintenance
- `check-status.sh` - Status check script for all services
- `README.md` - This documentation file

### GitHub Actions
- `.github/workflows/deploy.yml` - GitHub Actions workflow for CI/CD

## Architecture

```
Internet → Nginx (Port 80/443) → Django App (Port 8000)
                                      ↓
                                 PostgreSQL (Docker)
                                      ↓
                                   Redis (Docker)
```

## Security Features

- UFW firewall configured
- Fail2ban for intrusion prevention
- Services bound to localhost only
- Non-root user execution
- Automatic security updates
- Log rotation configured

## Manual Operations

### Check Status
```bash
./deploy/check-status.sh
```

### View Logs
```bash
# Application logs
sudo journalctl -u dicty-app -f

# Nginx logs
sudo tail -f /var/log/nginx/error.log
sudo tail -f /var/log/nginx/access.log

# Docker container logs
docker logs dicty-postgres -f
docker logs dicty-redis -f
```

### Restart Services
```bash
sudo systemctl restart dicty-app
sudo systemctl restart nginx
docker compose -f deploy/docker-compose.prod.yml restart
```

### Database Operations
```bash
# Connect to database
docker exec -it dicty-postgres psql -U resolwe -d resolwe

# Backup database
docker exec dicty-postgres pg_dump -U resolwe resolwe > backup_$(date +%Y%m%d).sql

# Restore database
cat backup.sql | docker exec -i dicty-postgres psql -U resolwe -d resolwe
```

### Manual Backup
```bash
sudo /usr/local/bin/dicty-backup.sh
```

## Troubleshooting

### Common Issues

1. **Service not starting**: Check logs with `journalctl -u dicty-app`
2. **Database connection issues**: Verify PostgreSQL container is running
3. **Frontend not loading**: Check nginx configuration and build files
4. **Permission issues**: Ensure dictyapp user owns project files

### Reset Deployment
```bash
# Stop all services
sudo systemctl stop dicty-app nginx
docker compose -f deploy/docker-compose.prod.yml down

# Clean up
docker system prune -f
sudo rm -rf /home/dictyapp/dicty_resolwe_server

# Re-run deployment
git clone https://github.com/ltrnavove/dicty_resolwe_server.git
cd dicty_resolwe_server
./deploy/deploy.sh
```

## Monitoring

The deployment includes:
- Daily automated backups
- Log rotation
- System resource monitoring
- Service health checks

## Environment Variables

Key environment variables (configured in `/home/dictyapp/dicty_resolwe_server/backend/.env`):
- `DEBUG=False` - Production mode
- `SECRET_KEY` - Django secret key (auto-generated)
- `ALLOWED_HOSTS` - Allowed hostnames
- `DB_*` - Database configuration
- `REDIS_*` - Redis configuration

## URLs

After deployment, the application will be available at:
- Frontend: http://95.179.242.134
- Django Admin: http://95.179.242.134/admin
- API: http://95.179.242.134/api/

## Support

For issues with deployment:
1. Check the status with `./deploy/check-status.sh`
2. Review logs with `journalctl -u dicty-app`
3. Verify all services are running
4. Check firewall settings with `ufw status`
