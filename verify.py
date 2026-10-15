#!/usr/bin/env python3
"""Independent exact checker. Standard library only; no search module imports.

The geometric first-crossing lemma is spelled out in REPORT.md. For every
pair of selected angular wedges, check an exact LP dual upper bound <=0
on simultaneous strict failure of the two radial faces. Arithmetic throughout
is fractions.Fraction, including square-root enclosures via integer isqrt.
"""
import argparse, gzip, hashlib, itertools, json, math
from fractions import Fraction as F
from pathlib import Path

if not __debug__:
    raise RuntimeError('Exact checks require Python without -O or -OO')

ZERO=F(0)

def dot(a,b): return sum((x*y for x,y in zip(a,b)),ZERO)
def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def det(a,b,c): return dot(cross(a,b),c)
def parse_cyl(c):
    assert len(c['vertices'])==3 and all(len(p)==3 for p in c['vertices'])
    assert len(c['axis'])==3
    return ([tuple(map(F,p)) for p in c['vertices']],tuple(map(F,c['axis'])))

def faces(c):
    """Return the two side faces and radial face as a.x<=b.
    Construct from 3D determinants, independently of angular formulas.
    """
    P,v=c; rows=[]
    for i,j,k in [(0,1,2),(2,0,1),(1,2,0)]:
        normal=cross(sub(P[j],P[i]),v); b=dot(normal,P[i])
        sign=dot(normal,P[k])-b
        assert sign!=0, 'Degenerate triangular prism'
        if sign>0: normal=tuple(-x for x in normal); b=-b
        scale=sum(map(abs,normal)); assert scale>0
        rows.append((tuple(x/scale for x in normal),b/scale))
    return rows

K=[((F(1),ZERO,F(1)),F(1)),((F(-1),ZERO,F(1)),F(1)),((ZERO,F(1),F(-1)),ZERO),((ZERO,F(-1),F(-1)),ZERO)]

def pair_rows(lf,uf):
    rows=[(a+(ZERO,),b) for a,b in K+lf[:2]+uf[:2]]
    for a,b in [lf[2],uf[2]]: rows.append((tuple(-x for x in a)+(F(1),),-b))
    return rows

def exact_solution(columns,target):
    """Solve an overdetermined full-column-rank system, checking consistency."""
    m=len(target); n=len(columns)
    assert 0<n<=m
    a=[[columns[j][i] for j in range(n)]+[target[i]] for i in range(m)]
    row=0; piv=[]
    for col in range(n):
        p=next((r for r in range(row,m) if a[r][col]),None)
        assert p is not None, 'Dependent dual support'
        a[row],a[p]=a[p],a[row]
        d=a[row][col]; a[row]=[x/d for x in a[row]]
        for r in range(m):
            if r!=row and a[r][col]:
                d=a[r][col]; a[r]=[x-d*y for x,y in zip(a[r],a[row])]
        piv.append(row); row+=1
    assert all(not any(a[r][:-1]) and a[r][-1]==0 for r in range(row,m)), 'Inconsistent exact dual equations'
    return [a[p][-1] for p in piv]

def check_dual(rows,basis,empty=False):
    assert len(set(basis))==len(basis) and all(0<=i<10 for i in basis)
    if empty:
        target=[ZERO]*4+[F(1)]
        cols=[list(rows[i][0])+[F(1)] for i in basis]
    else:
        target=[ZERO,ZERO,ZERO,F(1)]
        cols=[list(rows[i][0]) for i in basis]
    y=exact_solution(cols,target)
    assert all(t>=0 for t in y), 'Negative exact dual multiplier'
    bound=sum((t*rows[i][1] for i,t in zip(basis,y)),ZERO)
    assert bound<0 if empty else bound<=0, 'Dual fails to exclude simultaneous radial failure'
    return bound

def rotation(c):
    ps,v=c
    return ([(-y,x,1-t) for x,y,t in ps],(-v[1],v[0],-v[2]))

def same_cylinder(c,d):
    p,v=c; q,w=d
    return set(p)==set(q) and cross(v,w)==(0,0,0)

def check_geometry(data):
    assert data['format']=='angular-triangles-v1'
    L=list(map(parse_cyl,data['lower'])); U=list(map(parse_cyl,data['upper']))
    n=len(L); assert len(U)==n and n>0
    for family,cs in [('lower',L),('upper',U)]:
        previous=None
        for j,(P,v) in enumerate(cs):
            if family=='lower':
                assert P[0]==(0,0,0) and v[0]==1
                assert all(p[0]==0 for p in P) and P[1][2]>0 and P[2][2]>0
                qs=[p[1]/p[2] for p in P[1:]]
                shifts=[v[1]-q*v[2] for q in qs]
            else:
                assert P[0]==(0,0,1) and v[1]==1
                assert all(p[1]==0 for p in P) and P[1][2]<1 and P[2][2]<1
                qs=[p[0]/(1-p[2]) for p in P[1:]]
                shifts=[v[0]+q*v[2] for q in qs]
            assert -1<=qs[0]<qs[1]<=1
            if j==0: assert qs[0]==-1 and shifts[0]==0
            else: assert previous==(qs[0],shifts[0]), 'Angular supporting planes do not match'
            previous=(qs[1],shifts[1])
        assert previous==(1,0)
    for j in range(n):
        assert same_cylinder(rotation(L[j]),U[n-1-j]), 'Rotation symmetry fails'
        assert same_cylinder(rotation(U[j]),L[j]), 'Rotation symmetry fails'
    return L,U

def orbit(j,k,n): return {(j,k),(k,n-1-j),(n-1-j,n-1-k),(n-1-k,j)}

def area_bounds(cylinders,digits=30):
    D=10**digits; lo=ZERO; hi=ZERO; relative=ZERO
    for P,v in cylinders:
        coeff=abs(det(sub(P[1],P[0]),sub(P[2],P[0]),v))/2
        norm=v[0]**2+v[1]**2+2*v[2]**2
        assert coeff>0 and norm>0
        k=math.isqrt(D*D*norm.denominator//norm.numerator)
        assert F(k*k,D*D)*norm<=1 and F((k+1)**2,D*D)*norm>1
        lo+=coeff*F(k,D); hi+=coeff*F(k+1,D)
        relative+=coeff/(max(abs(v[0]),abs(v[2]))+max(abs(v[1]),abs(v[2])))
    return lo,hi,relative

def decimal_bound(x,digits=18,upper=False):
    z=x*10**digits; k=z.numerator//z.denominator
    if upper and F(k)!=z: k+=1
    return ('-' if k<0 else '')+str(abs(k)//10**digits)+'.'+str(abs(k)%10**digits).zfill(digits)

def check_euclidean(path,raw,L,U):
    data=json.loads(Path(path).read_text())
    assert data['format']=='rational-euclidean-triangular-cylinder-cover-v1'
    assert data['source_cover_sha256']==hashlib.sha256(raw).hexdigest()
    assert data['tetrahedron_vertices']==[[1,1,-1],[-1,-1,-1],[1,-1,1],[-1,1,1]]
    assert data['minimum_projection_area']=='2*sqrt(2)'
    assert len(data['cylinders'])==len(L)+len(U)
    for idx,((P,v),c) in enumerate(zip(L+U,data['cylinders'])):
        assert c['family']==('lower' if idx<len(L) else 'upper') and c['index']==idx%len(L)
        V=(v[0]+v[1],v[0]-v[1],2*v[2])
        Q=[(x+y,x-y,2*t-1) for x,y,t in P]
        assert tuple(map(F,c['axis']))==V
        assert [tuple(map(F,p)) for p in c['intercept_triangle']]==Q
        bases=[tuple(map(F,p)) for p in c['perpendicular_base']]
        assert len(bases)==3 and all(len(b)==3 for b in bases)
        for p,b in zip(Q,bases):
            assert dot(b,V)==0
            assert cross(sub(p,b),V)==(0,0,0)
        assert det(sub(Q[1],Q[0]),sub(Q[2],Q[0]),V)!=0
        normal=cross(sub(bases[1],bases[0]),sub(bases[2],bases[0]))
        coeff=abs(det(sub(P[1],P[0]),sub(P[2],P[0]),v))/2
        norm=v[0]**2+v[1]**2+2*v[2]**2
        # Euclidean base area squared / (2*sqrt(2))^2.
        assert dot(normal,normal)/32==coeff**2/norm
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def run(cover_path,cert_path,euclidean_path=None):
    raw=Path(cover_path).read_bytes(); data=json.loads(raw); L,U=check_geometry(data); n=len(L)
    lf=list(map(faces,L)); uf=list(map(faces,U))
    opener=gzip.open if str(cert_path).endswith('.gz') else open
    with opener(cert_path,'rt') as f: cert=json.load(f)
    assert cert['format']=='angular-pair-dual-supports-v1'
    assert cert['cover_sha256']==hashlib.sha256(raw).hexdigest(), 'Cover/certificate hash mismatch'
    expected={min(orbit(j,k,n)) for j in range(n) for k in range(n)}
    seen=set(); empty_count=0
    for record in cert['cells']:
        j,k,basis,empty=record
        assert (j,k) in expected and (j,k) not in seen
        check_dual(pair_rows(lf[j],uf[k]),basis,empty)
        seen.add((j,k)); empty_count+=int(empty)
    assert seen==expected, 'Incomplete angular-pair certificate'
    euclidean_hash=check_euclidean(euclidean_path,raw,L,U) if euclidean_path else None
    lo,hi,rel=area_bounds(L+U)
    assert hi<F(1,2), 'No certified saving below one half'
    result={'verified':True,'cylinders':2*n,'angular_pairs':n*n,'exact_dual_certificates':len(seen),'empty_cells':empty_count,
            'normalized_area_lower':decimal_bound(lo),'normalized_area_upper':decimal_bound(hi,upper=True),
            'saving_lower':decimal_bound(F(1,2)-hi),'directionwise_relative_area_upper':decimal_bound(rel,upper=True),
            'target_1e-3_reached':hi<=F(499,1000),'cover_sha256':hashlib.sha256(raw).hexdigest()}
    if euclidean_hash: result['rational_euclidean_cover_sha256']=euclidean_hash
    print(json.dumps(result,indent=2)); return result

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('cover',nargs='?',default='cover.json'); p.add_argument('certificate',nargs='?',default='certificate.json.gz'); p.add_argument('--output'); p.add_argument('--euclidean'); a=p.parse_args()
    result=run(a.cover,a.certificate,a.euclidean)
    if a.output: Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
