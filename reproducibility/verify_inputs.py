from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parent
manifest=json.loads((root/'manifest.json').read_text(encoding='utf-8'))
fail=[]
for row in manifest:
 p=root/row['path']
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:fail.append(row['path'])
print(json.dumps({'checked':len(manifest),'mismatches':fail},ensure_ascii=False))
raise SystemExit(bool(fail))
