# Publication status of the m6465 report payload

`m6465_report_payload.json` is the exact body POSTed to `$TENSOR_API_BASE/tensor/reports`
for mission 6465 (title, summary, six-element markdown, embedded figure). The endpoint
answered **HTTP 500** on every attempt (3 x, two paths, fresh admin token):

```
(psycopg.errors.UniqueViolation) duplicate key value violates unique constraint
"uq_reports_mission_md5"
DETAIL: Key (mission_id, md5(markdown))=(6458, c58298981ee9a89779d51247acff2ee9) already exists
```

The offending mission id is 6458, not 6465, and the markdown hash is not this payload's
(two different bodies produced the same key), so the server is writing a stale queued
report instead of the request body. That is a platform-side defect outside this mission's
scope (TensorPlatform is not editable by agents); the payload is committed here so the
watcher / 补报 mission can publish it verbatim once the endpoint is fixed:

```bash
curl -s -X POST $TENSOR_API_BASE/tensor/reports \
  -H "Authorization: Bearer $TOKEN" -H 'Content-Type: application/json' \
  -d @docs/intersect_resolution/m6465_report_payload.json
```
