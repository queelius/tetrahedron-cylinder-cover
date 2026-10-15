#!/usr/bin/env python3
"""Exact lower bound for symmetric even-polynomial secant angular ansatz.
See REPORT.md for the necessary diagonal-test-point argument. This does NOT
bound arbitrary cylinder covers or arbitrary direction profiles.
"""
import argparse,json,math,time
from fractions import Fraction as F
from pathlib import Path

if not __debug__:
    raise RuntimeError('Exact checks require Python without -O or -OO')

def prove(n=512,subdivisions=16,target=F(499,1000),shape=F(0),max_tilt=F(2)):
 assert n>0 and n%2==0 and subdivisions>0 and max_tilt>0
 coeffs=shape if isinstance(shape,list) else [shape]
 assert max_tilt<=1 or (all(c==0 for c in coeffs) and max_tilt<=2)
 D=10**20
 bins=[]
 for j in range(n//2,n):
  l=F(2*j,n)-1; r=l+F(2,n)
  pl=(1-l*l)/4*(1+sum(c*l**(2*i+2) for i,c in enumerate(coeffs)))
  pr=(1-r*r)/4*(1+sum(c*r**(2*i+2) for i,c in enumerate(coeffs)))
  assert pl>=0 and pr>=0
  a=(r*pl-l*pr)/(r-l); b=(pl-pr)/(r-l)
  assert 0<=b<1 and 0<a<1 and 0<=a-b<=1-r
  for k in range(subdivisions):
   q=l+(r-l)*F(k,subdivisions)
   # Tcrit is increasing in q>=0 and x=tau^2>=0.
   bins.append((a*a,a*b*q,a*a+2*b*b))
 width=F(2,n*subdivisions)
 def lower(L,U):
  s=0
  for aa,abq,A in bins:
   T=(1+L*aa)/(2*(1+L*(aa-abq)))
   z=1+U*A
   root=math.isqrt((D*D*z.numerator)//z.denominator)
   if F(root*root,D*D)<z: root+=1
   # root/D >= sqrt(z); round each contribution downward to grid 1/D.
   value=T*T*D*D/root
   s+=value.numerator//value.denominator
  return 2*width*F(s,D)
 pending=[(F(0),max_tilt**2,0)]; accepted=[]; count=0; start=time.time()
 while pending:
  L,U,depth=pending.pop(); bound=lower(L,U); count+=1
  if bound>=target:
   accepted.append([str(L),str(U),str(bound)])
  else:
   assert depth<20, f'Cannot prove requested bound on [{L},{U}]: {float(bound)}'
   M=(L+U)/2; pending.extend([(M,U,depth+1),(L,M,depth+1)])
  if count%100==0: print(f'intervals tested {count}, accepted {len(accepted)}, pending {len(pending)}, elapsed {time.time()-start:.1f}s',flush=True)
 result={'verified':True,'scope':'same even radial profile in both families; n equal secant angular sectors',
         'shape':[str(c) for c in coeffs],'max_tilt':str(max_tilt),
         'n':n,'q_subdivisions_per_sector':subdivisions,'normalized_area_lower':str(target),
         'intervals_tested':count,'accepted_intervals':accepted,'seconds':time.time()-start}
 print(json.dumps({k:v for k,v in result.items() if k!='accepted_intervals'},indent=2))
 return result
if __name__=='__main__':
 p=argparse.ArgumentParser(); p.add_argument('--n',type=int,default=512); p.add_argument('--subdivisions',type=int,default=16); p.add_argument('--target',default='499/1000'); p.add_argument('--out',default='structural_bound.json');p.add_argument('--shape',default='0');p.add_argument('--max-tilt',default='2');a=p.parse_args()
 Path(a.out).write_text(json.dumps(prove(a.n,a.subdivisions,F(a.target),[F(x) for x in a.shape.split(',')],F(a.max_tilt)),indent=2)+'\n')
