"""Independent exact angular cell decomposition and primal polyhedron clipping.
Does not read the supplied dual certificate. Each pair cell is clipped by both
radial failure halfspaces; absence of full dimension proves no open gap.
"""
from independent import *
import multiprocessing as mp,os

def proportional(a,b):
 k=next(i for i,z in enumerate(a) if z)
 r=b[k]/a[k]
 assert all(y==r*x for x,y in zip(a,b))
 return r

def check_fan(hs):
 body=initial()[0]
 assert all(ev(hs[0][0],p)<=0 for p in body)
 assert all(ev(hs[-1][2],p)<=0 for p in body)
 for h,g in zip(hs,hs[1:]):assert proportional(h[2],g[0])<0
 # If s_j is positive on the inside of the first boundary, we have s_0>=0,
 # s_n<=0. A first transition from >=0 to <=0 selects both angular sides.

def rem(hs):
 p=initial()
 for h in [hs[0],hs[2],tuple(-z for z in hs[1])]:
  p=clip(p,h)
  if p is None:break
 return p

UP=None;LOW=None

def row(i):
 p=LOW[i]; counts=[0,0,0]; maxbits=0
 if p is None:return i,counts,None
 # A pair cell which fails both cylinders has the following closed polytope.
 # Clip rather than use any LP solver, support, multiplier, or dual certificate.
 for j,hs in enumerate(UP):
  q=p
  for stage,h in enumerate([hs[0],hs[2],tuple(-z for z in hs[1])]):
   q=clip(q,h)
   if q is None:counts[stage]+=1;break
  if q is not None:
   w=mul(tuple(map(sum,zip(*q[0]))),Q(1,len(q[0])))
   return i,counts,{'pair':[i,j],'point':list(map(str,w))}
 return i,counts,None

def main():
 global UP,LOW
 ap=argparse.ArgumentParser();ap.add_argument('--rows',type=int,default=2048);ap.add_argument('--workers',type=int,default=3);args=ap.parse_args()
 d=json.load(open('../cover.json'))
 hsL=[planes(c) for c in d['lower']];UP=[planes(c) for c in d['upper']]
 check_fan(hsL);check_fan(UP)
 LOW=[rem(hs) for hs in hsL]
 t=time.monotonic();counts=[0,0,0];done=0
 # fork shares read-only exact geometry; only three workers perform computation.
 with mp.get_context('fork').Pool(args.workers) as pool:
  for i,c,bad in pool.imap_unordered(row,range(args.rows),chunksize=1):
   counts=[a+b for a,b in zip(counts,c)];done+=1
   if bad:
    open('gap.json','w').write(json.dumps(bad,indent=2)+'\n');raise RuntimeError(str(bad))
   if done%8==0 or done==args.rows:
    r={'rows_done':done,'pair_cells':sum(counts),'empty_after_clip':counts,'seconds':round(time.monotonic()-t,3),'pid':os.getpid()}
    print(json.dumps(r),flush=True)
    open('cells-progress.json','w').write(json.dumps(r,indent=2)+'\n')
 r['all_cells_checked']=done==2048 and sum(counts)==2048**2
 open('cells-result.json','w').write(json.dumps(r,indent=2)+'\n')
if __name__=='__main__':main()
