cd /Users/lenatrnovec/backend/dicty_resolwe_server/resolwe_server && python manage.py runworker resolwe-server.manager.control

pkill -f runworker || true
pkill -f runlistener || true
cd /Users/lenatrnovec/backend/dicty_resolwe_server/resolwe_server && python manage.py runlistener

cd /Users/lenatrnovec/backend/dicty_resolwe_server && python setup_annotations.py