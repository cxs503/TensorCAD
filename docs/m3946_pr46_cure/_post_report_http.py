import json, urllib.request
from pathlib import Path

token = json.load(open('/tmp/login_out.json'))['access_token']
body = Path('/tmp/report.json').read_bytes()
base = 'http://127.0.0.1:8010/api/v1'
req = urllib.request.Request(base + '/tensor/reports', data=body, method='POST',
                             headers={'Authorization': 'Bearer ' + token,
                                      'Content-Type': 'application/json'})
try:
    with urllib.request.urlopen(req, timeout=120) as r:
        print('HTTP', r.status, r.read()[:300])
except urllib.error.HTTPError as e:
    print('HTTPError', e.code, e.read()[:300])