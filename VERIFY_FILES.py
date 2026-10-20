from pathlib import Path
import csv, hashlib
root=Path(__file__).resolve().parent
with (root/'SHA256SUMS.csv').open(encoding='utf-8-sig',newline='') as f:
    rows=list(csv.DictReader(f))
for row in rows:
    p=root/row['file']
    assert p.is_file(),f"Missing {row['file']}"
    assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'],f"Changed {row['file']}"
print(f'PASS: {len(rows)} distributed files match SHA256SUMS.csv')
