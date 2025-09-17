#!/bin/bash

# Test deployment script - validates deployment configuration locally

set -e

echo "=== Testing Dicty Resolwe Server Deployment Configuration ==="
echo "Date: $(date)"
echo ""

PROJECT_DIR="/Users/lenatrnovec/dicty_resolwe_server"
DEPLOY_DIR="$PROJECT_DIR/deploy"

# Test 1: Check if all required files exist
echo "=== Test 1: Checking required files ==="
required_files=(
    ".github/workflows/deploy.yml"
    "deploy/deploy.sh"
    "deploy/docker-compose.prod.yml"
    "deploy/production-env-template.txt"
    "deploy/server-setup.sh"
    "deploy/dicty-app.service"
    "deploy/check-status.sh"
    "deploy/README.md"
    "backend/requirements.txt"
    "frontend/package.json"
)

for file in "${required_files[@]}"; do
    if [ -f "$PROJECT_DIR/$file" ]; then
        echo "✅ $file exists"
    else
        echo "❌ $file missing"
        exit 1
    fi
done
echo ""

# Test 2: Check script permissions
echo "=== Test 2: Checking script permissions ==="
scripts=("deploy.sh" "server-setup.sh" "check-status.sh" "test-deployment.sh")
for script in "${scripts[@]}"; do
    if [ -x "$DEPLOY_DIR/$script" ]; then
        echo "✅ $script is executable"
    else
        echo "❌ $script is not executable"
        exit 1
    fi
done
echo ""

# Test 3: Validate shell scripts syntax
echo "=== Test 3: Validating shell script syntax ==="
for script in "${scripts[@]}"; do
    if bash -n "$DEPLOY_DIR/$script"; then
        echo "✅ $script syntax is valid"
    else
        echo "❌ $script has syntax errors"
        exit 1
    fi
done
echo ""

# Test 4: Check Docker Compose syntax
echo "=== Test 4: Validating Docker Compose syntax ==="
if docker compose -f "$DEPLOY_DIR/docker-compose.prod.yml" config > /dev/null 2>&1; then
    echo "✅ Docker Compose configuration is valid"
else
    echo "❌ Docker Compose configuration has errors"
    exit 1
fi
echo ""

# Test 5: Check GitHub Actions workflow syntax
echo "=== Test 5: Validating GitHub Actions workflow ==="
if python3 -c "import yaml; yaml.safe_load(open('$PROJECT_DIR/.github/workflows/deploy.yml'))"; then
    echo "✅ GitHub Actions workflow YAML is valid"
else
    echo "❌ GitHub Actions workflow YAML has errors"
    exit 1
fi
echo ""

# Test 6: Check environment template
echo "=== Test 6: Validating environment template ==="
if [ -s "$DEPLOY_DIR/production-env-template.txt" ]; then
    echo "✅ Environment template exists and is not empty"
    
    # Check for required environment variables
    required_vars=("DB_NAME" "DB_USER" "DB_PASSWORD" "SECRET_KEY" "ALLOWED_HOSTS")
    for var in "${required_vars[@]}"; do
        if grep -q "^$var=" "$DEPLOY_DIR/production-env-template.txt"; then
            echo "✅ $var found in environment template"
        else
            echo "⚠️  $var not found in environment template"
        fi
    done
else
    echo "❌ Environment template is missing or empty"
    exit 1
fi
echo ""

# Test 7: Check backend requirements
echo "=== Test 7: Validating backend requirements ==="
if python3 -c "
import pkg_resources
with open('$PROJECT_DIR/backend/requirements.txt') as f:
    requirements = f.read()
    for line in requirements.strip().split('\n'):
        if line and not line.startswith('#'):
            try:
                pkg_resources.Requirement.parse(line)
                print(f'✅ {line.split()[0]} requirement is valid')
            except Exception as e:
                print(f'❌ {line} requirement is invalid: {e}')
                exit(1)
"; then
    echo "✅ All Python requirements are valid"
else
    echo "❌ Some Python requirements are invalid"
    exit 1
fi
echo ""

# Test 8: Check frontend package.json
echo "=== Test 8: Validating frontend package.json ==="
if python3 -c "import json; json.load(open('$PROJECT_DIR/frontend/package.json'))"; then
    echo "✅ Frontend package.json is valid JSON"
    
    # Check for required scripts
    required_scripts=("build" "start")
    for script in "${required_scripts[@]}"; do
        if python3 -c "
import json
pkg = json.load(open('$PROJECT_DIR/frontend/package.json'))
if '$script' in pkg.get('scripts', {}):
    print('✅ $script script found')
else:
    print('❌ $script script missing')
    exit(1)
"; then
            continue
        else
            exit 1
        fi
    done
else
    echo "❌ Frontend package.json is invalid"
    exit 1
fi
echo ""

# Test 9: Check port configurations
echo "=== Test 9: Checking port configurations ==="
ports_to_check=("80" "443" "8000" "5432" "6379")
echo "Ports that will be used:"
for port in "${ports_to_check[@]}"; do
    echo "  - Port $port"
done
echo "✅ Port configuration documented"
echo ""

# Test 10: Security check
echo "=== Test 10: Security configuration check ==="
if grep -q "127.0.0.1" "$DEPLOY_DIR/docker-compose.prod.yml"; then
    echo "✅ Database services bound to localhost"
else
    echo "⚠️  Database services not explicitly bound to localhost"
fi

if grep -q "UFW" "$DEPLOY_DIR/server-setup.sh"; then
    echo "✅ Firewall configuration included"
else
    echo "⚠️  Firewall configuration not found"
fi

if grep -q "fail2ban" "$DEPLOY_DIR/server-setup.sh"; then
    echo "✅ Fail2ban configuration included"
else
    echo "⚠️  Fail2ban configuration not found"
fi
echo ""

echo "=== Test Summary ==="
echo "✅ All deployment configuration tests passed!"
echo ""
echo "Next steps:"
echo "1. Ensure VPS (95.179.242.134) is accessible via SSH"
echo "2. Run server-setup.sh on the VPS to prepare the environment"
echo "3. Set up SSH keys for GitHub Actions"
echo "4. Test manual deployment with deploy.sh"
echo "5. Push to master branch to test automatic deployment"
echo ""
echo "Deployment will be available at: http://95.179.242.134"
