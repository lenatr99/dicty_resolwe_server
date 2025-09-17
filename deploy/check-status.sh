#!/bin/bash

# Status check script for Dicty Resolwe Server

echo "=== Dicty Resolwe Server Status Check ==="
echo "Date: $(date)"
echo ""

# Check system resources
echo "=== System Resources ==="
echo "Memory usage:"
free -h
echo ""
echo "Disk usage:"
df -h /
echo ""
echo "Load average:"
uptime
echo ""

# Check Docker containers
echo "=== Docker Containers ==="
docker ps --format "table {{.Names}}\t{{.Status}}\t{{.Ports}}"
echo ""

# Check services
echo "=== System Services ==="
services=("dicty-app" "nginx" "docker" "fail2ban")
for service in "${services[@]}"; do
    if systemctl is-active --quiet $service; then
        echo "✅ $service: running"
    else
        echo "❌ $service: not running"
    fi
done
echo ""

# Check network connectivity
echo "=== Network Connectivity ==="
if curl -f http://localhost:8000/api/ > /dev/null 2>&1; then
    echo "✅ Backend API: responding"
else
    echo "❌ Backend API: not responding"
fi

if curl -f http://localhost/ > /dev/null 2>&1; then
    echo "✅ Frontend: responding"
else
    echo "❌ Frontend: not responding"
fi

# Check database connectivity
echo ""
echo "=== Database Connectivity ==="
if docker exec dicty-postgres pg_isready -U resolwe > /dev/null 2>&1; then
    echo "✅ PostgreSQL: ready"
else
    echo "❌ PostgreSQL: not ready"
fi

if docker exec dicty-redis redis-cli ping > /dev/null 2>&1; then
    echo "✅ Redis: responding"
else
    echo "❌ Redis: not responding"
fi

# Check logs for errors
echo ""
echo "=== Recent Errors ==="
echo "Last 5 error lines from dicty-app service:"
journalctl -u dicty-app --no-pager -n 5 -p err --since "1 hour ago" || echo "No recent errors"

echo ""
echo "Last 5 error lines from nginx:"
tail -n 5 /var/log/nginx/error.log 2>/dev/null || echo "No nginx error log found"

echo ""
echo "=== Firewall Status ==="
ufw status

echo ""
echo "=== Port Status ==="
echo "Listening ports:"
netstat -tlnp | grep -E ":(80|443|8000|5432|6379)" || echo "No relevant ports found"

echo ""
echo "=== Application URLs ==="
echo "Frontend: http://95.179.242.134"
echo "Django Admin: http://95.179.242.134/admin"
echo "API: http://95.179.242.134/api/"
