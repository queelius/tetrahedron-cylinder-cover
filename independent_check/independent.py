"""Independent exact geometric subtraction; no candidate certificate or code used.
All geometric predicates use FLINT rationals. Plane convention: a.p+d <= 0.
"""
import json,sys,time,hashlib,argparse
from itertools import combinations
from flint import fmpq as Q, arb, ctx

def vec(v): return tuple(Q(str(x)) for x in v)
def add(a,b): return tuple(x+y for x,y in zip(a,b))
def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def mul(a,s): return tuple(x*s for x in a)
def dot(a,b): return sum((x*y for x,y in zip(a,b)),Q(0))
def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def ev(h,v): return dot(h[:3],v)+h[3]
def planes(c):
 v=list(map(vec,c['vertices'])); a=vec(c['axis']); hs=[]
 for i,j,k in [(0,1,2),(1,2,0),(2,0,1)]:
  n=cross(sub(v[j],v[i]),a); h=(*n,-dot(n,v[i]))
  if ev(h,v[k])>0: h=tuple(-z for z in h)
  assert ev(h,v[k])<0
  hs.append(h)
 return hs

# Convex polyhedron: extreme vertices and unordered vertex sets of facets.
# Edges are intersections of distinct facets. All clipping is exact.
def dimension3(v):
 if len(v)<4:return False
 a=v[0]; b=next((sub(w,a) for w in v[1:] if w!=a),None)
 if b is None:return False
 n=next((cross(b,sub(w,a)) for w in v[1:] if any(cross(b,sub(w,a)))),None)
 return n is not None and any(dot(n,sub(w,a)) for w in v[1:])
def clip(poly,h,vals=None):
 v,fs=poly
 s=vals if vals is not None else [ev(h,p) for p in v]
 if all(z<=0 for z in s):return poly
 if all(z>=0 for z in s):return None # retain only full-dimensional pieces
 edges=set()
 for f,g in combinations(fs,2):
  e=f&g
  if len(e)==2: edges.add(tuple(sorted(e)))
  elif len(e)>2: raise AssertionError('non-extreme edge vertices')
 nv=[]; old={}; cut={}
 for i,z in enumerate(s):
  if z<=0:old[i]=len(nv);nv.append(v[i])
 for i,j in edges:
  if s[i]*s[j]<0:
   p=add(v[i],mul(sub(v[j],v[i]),s[i]/(s[i]-s[j])))
   cut[i,j]=len(nv);nv.append(p)
 if not dimension3(nv):return None
 nfs=[]
 for f in fs:
  nf={old[i] for i in f if i in old}
  nf.update(idx for (i,j),idx in cut.items() if i in f and j in f)
  if len(nf)>=3:nfs.append(frozenset(nf))
 cap={old[i] for i,z in enumerate(s) if z==0}
 cap.update(cut.values())
 if len(cap)>=3:nfs.append(frozenset(cap))
 assert len(nfs)>=4
 return nv,list(dict.fromkeys(nfs))

def subtract(poly,hs):
 # Fast rejection avoids clipping when a single facet separates the cylinder.
 values=[[ev(h,p) for p in poly[0]] for h in hs]
 if any(all(z>=0 for z in row) for row in values):return [poly]
 if all(all(z<=0 for z in row) for row in values):return []
 out=[];inside=poly
 for h in hs:
  s=[ev(h,p) for p in inside[0]]
  if all(z<=0 for z in s):continue
  if all(z>=0 for z in s):out.append(inside);return out
  outside=clip(inside,tuple(-z for z in h),[-z for z in s])
  if outside:out.append(outside)
  inside=clip(inside,h,s)
  if inside is None:return out
 return out

def initial():
 v=list(map(vec,[(1,0,0),(-1,0,0),(0,1,1),(0,-1,1)]))
 return v,[frozenset(f) for f in combinations(range(4),3)]

def geometry_checks(d,e):
 # Similarity F(x,y,t)=(x+y,x-y,2t-1), factor sqrt(2).
 def F(v): x,y,t=v;return (x+y,x-y,2*t-1)
 def D(v): x,y,t=v;return (x+y,x-y,2*t)
 tv=list(map(vec,e['tetrahedron_vertices']))
 assert set(tv)==set(map(F,initial()[0]))
 assert {dot(sub(a,b),sub(a,b)) for a,b in combinations(tv,2)}=={Q(8)}
 cs=d['lower']+d['upper'];assert len(cs)==len(e['cylinders'])==4096
 ctx.prec=192; total=arb(0); areas=[]
 for i,(c,z) in enumerate(zip(cs,e['cylinders'])):
  a=D(vec(c['axis'])); v=list(map(F,map(vec,c['vertices'])));b=list(map(vec,z['perpendicular_base']))
  assert a==vec(z['axis']) and v==list(map(vec,z['intercept_triangle']))
  assert a!=(0,0,0)
  for w,p in zip(v,b):
   assert dot(a,p)==0
   assert p==sub(w,mul(a,dot(a,w)/dot(a,a)))
  # Exact squared area from perpendicular base, and independently projected intercept.
  cp=cross(sub(b[1],b[0]),sub(b[2],b[0]));sq=dot(cp,cp)/4
  n=cross(sub(v[1],v[0]),sub(v[2],v[0]));sq2=dot(n,a)**2/(4*dot(a,a))
  assert sq==sq2 and sq>0
  areas.append(str(sq))
  total+=arb(sq).sqrt()
 cost=total/(2*arb(2).sqrt())
 assert cost.upper()<arb('0.499613830783714149').lower()
 result={'count':len(cs),'precision_bits':ctx.prec,'normalized_area_ball':str(cost),'total_base_area_ball':str(total),'saving_ball':str(arb('0.5')-cost),'area_squares_sha256':hashlib.sha256('\n'.join(areas).encode()).hexdigest()}
 open('area_squares.json','w').write(json.dumps(areas,separators=(',',':'))+'\n')
 return result

def coverage(cs,limit=None):
 polys=[initial()];t=time.monotonic();maxpieces=1
 for i,c in enumerate(cs[:limit]):
  hs=planes(c);out=[]
  for p in polys:out.extend(subtract(p,hs))
  polys=out;maxpieces=max(maxpieces,len(polys))
  if i%32==0 or not polys:
   print(json.dumps({'step':i+1,'pieces':len(polys),'maxpieces':maxpieces,'seconds':round(time.monotonic()-t,2)}),flush=True)
  if not polys:break
 return polys,{'processed':i+1,'remaining_pieces':len(polys),'maxpieces':maxpieces,'seconds':time.monotonic()-t}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--limit',type=int);ap.add_argument('--areas-only',action='store_true');args=ap.parse_args()
 d=json.load(open('../cover.json'));e=json.load(open('../cover_euclidean.json'))
 result={'areas':geometry_checks(d,e)};print(json.dumps(result),flush=True)
 if not args.areas_only:
  rest,result['coverage']=coverage(d['lower']+d['upper'],args.limit)
  if rest:
   w=mul(tuple(map(sum,zip(*rest[0][0]))),Q(1,len(rest[0][0])))
   result['remainder_centroid']=list(map(str,w))
  open('independent-result.json','w').write(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
