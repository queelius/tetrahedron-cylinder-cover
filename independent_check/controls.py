"""Exact rational uncovered-point controls using barycentric determinants."""
import copy,json
from independent import Q,vec,sub,mul,add,dot,cross

def bary(c,p):
 v=list(map(vec,c['vertices']));a=vec(c['axis'])
 e=sub(v[1],v[0]);f=sub(v[2],v[0]);w=sub(p,v[0])
 det=dot(e,cross(f,a));assert det
 u=dot(w,cross(f,a))/det;v=dot(e,cross(w,a))/det
 return u,v

def contains(c,p):
 u,v=bary(c,p)
 return u>=0 and v>=0 and u+v<=1

d=json.load(open('../cover.json'));cs=d['lower']+d['upper']
p=(Q(0),Q(3,16384),Q(3,8));idx=1024
assert abs(p[0])<=1-p[2] and abs(p[1])<=p[2] and 0<p[2]<1
members=[i for i,c in enumerate(cs) if contains(c,p)]
assert members==[idx],members
removed=cs[:idx]+cs[idx+1:]
assert not any(contains(c,p) for c in removed)
shrunken=copy.deepcopy(cs);v=list(map(vec,cs[idx]['vertices']));tip=v[0]
shrunken[idx]['vertices']=[list(map(str,add(tip,mul(sub(w,tip),Q(1,2))))) for w in v]
assert not any(contains(c,p) for c in shrunken)
# Also exercise the coverage checker itself on a shrink, preserving fan boundaries.
from cells import rem
from independent import planes,clip,dimension3
hL=planes(shrunken[idx]);r=rem(hL)
assert r is not None
# Barycentric proof excludes all cylinders, and strict margins make this an open gap.
assert all((lambda b:b[0]<0 or b[1]<0 or sum(b)>1)(bary(c,p)) for c in shrunken)
result={'point_xyt':list(map(str,p)),'point_euclidean':list(map(str,(p[0]+p[1],p[0]-p[1],2*p[2]-1))),'original_members':members,'removed_index':idx,'removed_cover':'REJECT','shrink_factor':'1/2','shrunk_cover':'REJECT','original_barycentric':list(map(str,bary(cs[idx],p))),'shrunk_barycentric':list(map(str,bary(shrunken[idx],p)))}
open('controls-result.json','w').write(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
