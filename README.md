cd /Users/lenatrnovec/dicty_resolwe_server/resolwe_server && python manage.py runworker resolwe-server.manager.control

pkill -f runworker || true
pkill -f runlistener || true
cd /Users/lenatrnovec/dicty_resolwe_server/resolwe_server && python manage.py runlistener

python setup_annotations.py