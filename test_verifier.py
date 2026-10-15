"""Adversarial controls for the exact checker, using only the standard library."""
import contextlib,copy,gzip,hashlib,io,json,tempfile
from fractions import Fraction as F
from pathlib import Path
import verify

raw=Path('cover.json').read_bytes(); data=json.loads(raw)
with gzip.open('certificate.json.gz','rt') as f: cert=json.load(f)

def rejected(d,c):
 with tempfile.TemporaryDirectory(dir='.') as td:
  root=Path(td); r=(json.dumps(d,separators=(',',':'))+'\n').encode(); (root/'cover.json').write_bytes(r)
  c=copy.deepcopy(c); c['cover_sha256']=hashlib.sha256(r).hexdigest()
  with gzip.open(root/'cert.json.gz','wt') as f: json.dump(c,f)
  try:
   with contextlib.redirect_stdout(io.StringIO()): verify.run(root/'cover.json',root/'cert.json.gz')
  except AssertionError: return
  raise AssertionError('Corrupted input was accepted')

# Exact area sanity check: two Bang triangles each contribute 1/4.
P=[(F(0),F(0),F(0)),(F(0),F(-1,2),F(1,2)),(F(0),F(1,2),F(1,2))]
lo,hi,rel=verify.area_bounds([(P,(F(1),F(0),F(0)))]*2)
assert lo<=F(1,2)<=hi and hi-lo<F(1,10**28) and rel==F(1,2)

bad=copy.deepcopy(data); bad['lower'][0]['axis'][1]=str(F(bad['lower'][0]['axis'][1])+F(1,100))
rejected(bad,cert)
badcert=copy.deepcopy(cert); badcert['cells'].pop(); rejected(data,badcert)
badcert=copy.deepcopy(cert); badcert['cells'][0][2]=[]; rejected(data,badcert)

# Preserve angular matching and rotation symmetry, but pull every radial face inward.
bad=copy.deepcopy(data)
for name in ['lower','upper']:
 for c in bad[name]:
  tip=list(map(F,c['vertices'][0]))
  for i in [1,2]: c['vertices'][i]=[str(o+(F(x)-o)*F(999,1000)) for x,o in zip(c['vertices'][i],tip)]
rejected(bad,cert)
euc=json.loads(Path('cover_euclidean.json').read_text())
euc['cylinders'][0]['perpendicular_base'][0][0]=str(F(euc['cylinders'][0]['perpendicular_base'][0][0])+1)
with tempfile.TemporaryDirectory(dir='.') as td:
 path=Path(td)/'bad_euclidean.json'; path.write_text(json.dumps(euc))
 L,U=verify.check_geometry(data)
 try: verify.check_euclidean(path,raw,L,U)
 except AssertionError: pass
 else: raise AssertionError('Corrupted perpendicular base was accepted')
print('PASS: exact Bang area control; broken angular matching, missing cell, invalid dual, inward radial perturbation, and corrupted Euclidean base all rejected.')
