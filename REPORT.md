# The structural barrier for the final symmetric family

Alexander Richard Towell, with substantial OpenAI Codex assistance.
This supplies the geometric argument referenced by structural_bound.py for
the final family at source commit a73f6955cdd97b4184081e03581bb6cb6190c33b.
Ricardo Pinho's construction, normalization and original arithmetic certificate
are the predecessors. This is a written proof with exact finite arithmetic,
not a formal proof or human peer review.

## Statement and restrictions

Let n=2048, q_j=-1+2j/n, and
phi(q)=(1-q^2)(1-q^2/200+27q^4/1250)/4.
On a sector [l,r], write its secant as phi_j(q)=a-bq, where
a=(r phi(l)-l phi(r))/(r-l), b=(phi(l)-phi(r))/(r-l).

Use precisely the lower and upper cylinders in README.md, with respective
axes (1,tau a,tau b) and (-tau a,1,tau b). The positive node radii R_j are
arbitrary subject to R_j=R_{n-j} and use of the same radii in both families.
For 0<=tau<=1, if these cylinders cover
K={|x|<=1-t, |y|<=t}, then their normalized total perpendicular-base area
satisfies C>=4993/10000.

This statement does not cover different axes, nonsymmetric radii, different
partitions, or all parameter choices accepted by structural_bound.py.
In particular, this note does not certify its default max_tilt=2 setting.
The explicit 4096-cylinder cover and its full pair-cell certificate are
separate results, not assumptions used here.

The finite geometric conditions needed below, for every positive-half
sector 0<=l<r<=1, are
  0<=phi(l), phi(r),  0<a<1,  0<=b<1,  0<=a-b<=1-r.
They are verified by rational arithmetic for these 1024 sectors.
Evenness supplies a unchanged and b negated on the reflected sector.
These conditions also imply 0<=phi(q_j)<1 at every node.

## 1. Radial interpolation and area

Put z_j=1/R_j and interpolate z linearly on each sector. The intercept
triangle conv{0,(l R_l,R_l),(r R_r,R_r)} in (y,t) coordinates contains
(q T,T), for q in [l,r] and T>=0, exactly when T z(q)<=1.
Indeed its two nonzero barycentric weights are
  T(r-q) z_l/(r-l),  T(q-l) z_r/(r-l).
Thus the ray cutoff is R(q)=1/z(q); it is not the linear interpolant of R_j.
The upper family has the same cutoff, and z(-q)=z(q).

The physical metric is diag(1,1,2). To check the normalization, the four
physical face area-vectors are (plus/minus sqrt(2),0,-1) and
(0,plus/minus sqrt(2),1). Half the sum of their absolute dot products
with a physical unit direction u gives projection area
  max(sqrt(2)|u_x|,|u_z|)+max(sqrt(2)|u_y|,|u_z|).
Set U=max(|u_x|,|u_z|/sqrt(2)), V=max(|u_y|,|u_z|/sqrt(2)). Then
(U+V)^2>=u_x^2+u_y^2+u_z^2=1, so this area is at least sqrt(2),
with equality along the x axis. Thus the minimum projection area is sqrt(2).
For one lower intercept triangle its cross product
is along the x axis, with magnitude sqrt(2)(r-l) R_l R_r. Its physical axis
has norm sqrt(1+tau^2(a^2+2b^2)); projection to the perpendicular plane
therefore gives normalized area
  (r-l) R_l R_r / (2 sqrt(1+tau^2(a^2+2b^2))).
The upper triangle has the same area. Hence
  C = sum_j (r-l)/(z_l z_r sqrt(1+tau^2(a^2+2b^2)))
    = sum_j integral_l^r R(q)^2/sqrt(1+tau^2(a^2+2b^2)) dq.
The integral identity follows directly by integrating the reciprocal square
of a positive affine function; it remains true when z_l=z_r.

## 2. A point that tests both reflected radial faces

Fix a positive-half sector and l<q<r, with q>0. Set s=tau^2,
  D = 1+s a(a-bq),    T = (1+s a^2)/(2D),
  x = -qT(1+tau a)/(1+s a^2),
  y =  qT(1-tau a)/(1+s a^2),
  t = T+tau b x.
Here D>=1, T>0, x<0, y>0, and |x|>=y because tau a<1.
Direct substitution gives
  y-tau a x=qT,       t-tau b x=T,
  x+tau a y=-qT,      1-t-tau b y=T.
The first pair is the lower intercept in the chosen sector. The second is
the upper intercept in the reflected sector, whose coefficients are a,-b.
Thus both intercepts have the same radial coordinate T and opposite q.

This point lies in K. Specifically,
  2D(t-y) = 1-q+tau q(a-b)+s a(a-bq) >= 0,
  2D(1-t+x) = 1-q-tau q(a-b)+s a(a-bq) >= 0.
For the second inequality use
tau q(a-b)<=a-b<=1-r<=1-q. Since x<0<y, these two inequalities
are exactly |x|<=1-t and |y|<=t.

## 3. Other angular sectors cannot supply coverage

This step is necessary: the two chosen radial faces alone would not suffice
if a different cylinder could cover the test point.

For lower cylinders define their common angular boundary heights
  H_k=q_k t+tau phi(q_k)x.
Membership in sector k requires H_k<=y<=H_{k+1}.
For nonpositive nodes H_k<=0<y. On positive-half sectors,
  H_{k+1}-H_k=(q_{k+1}-q_k)(t-tau b_k x)>0,
since t>=y>0, x<0 and b_k>=0.
For the chosen sector, substitution gives
  y-H_j=(q-l)T>0,    H_{j+1}-y=(r-q)T>0.
It is therefore the only possible lower sector.

For upper cylinders define
  G_k=q_k(1-t)-tau phi(q_k)y.
Membership requires G_k<=x<=G_{k+1}.
At every nonnegative node,
  G_k>=-tau phi(q_k)y>-y>=x,
using 0<=phi(q_k)<1 and tau<=1.
On negative-half sectors b_k is negative and |b_k|<1, so
  G_{k+1}-G_k=(q_{k+1}-q_k)(1-t+tau b_k y)>0,
because 1-t>=|x|>=y>0.
The reflected-sector intercept equations place x strictly between its
two G boundaries. It is the only possible upper sector.

Consequently full coverage forces the lower sector or its upper reflection
to contain the point. Both conditions are identical, namely T z(q)<=1.
Thus, on every positive-half sector,
  R(q) >= T(q,s) = (1+s a^2)/(2(1+s(a^2-abq))).

This proves the inequality for sector interiors. Positivity and continuity
of z give the one-sided limits at all sector endpoints (including q=0,1).
Even where the two sector formulas have different limiting values, each
one-sided constraint follows separately. These finitely many endpoints
also have no effect on the integral. The derivation permits tau=0 and
tau=1 throughout; neither is handled by an unproved parameter limit.

## 4. The integral bound and exact arithmetic

Reflection in q changes b to -b and preserves a,z. Therefore
  C >= 2 sum_{j=n/2}^{n-1} integral_l^r
       T(q,s)^2 / sqrt(1+s(a^2+2b^2)) dq.

The relevant derivatives have nonnegative numerators:
  dT/dq = s ab(1+s a^2)/(2D^2),
  dT/ds = abq/(2D^2).
On a squared-tilt interval L<=s<=U and a q subinterval [u,v] in the
chosen sector, its integral is consequently at least
  (v-u) T(u,L)^2 / sqrt(1+U(a^2+2b^2)).

The original script uses four q subintervals per sector, upper bounds
each square root, rounds the contribution down, and bisects the
squared-tilt domain. Its saved structural_final_family.json reports
210 accepted intervals, from 419 binary-tree visits, covering [0,1].

The independent checker supplied with this note verifies the complete
interval partition and fixed-family metadata, checks all geometric
coefficient inequalities, and recomputes a lower enclosure on every
interval. It does not use the JSON's claimed lower-bound values as proof.
Its arithmetic procedure and actual verification results are documented
with the checker. Passing those checks certifies the last finite step
of this written argument; running the old script alone does not supply
Sections 1–3.

Run from the repository root (Python standard library only):

```sh
python structural_independent.py --output structural_independent_result.json
python -m unittest -v test_structural_independent.py
```

The independent run checks all 1024 coefficient pairs and all 860160
q/squared-tilt cells in the 210 intervals. Its minimum lower enclosure is
71956848089952005/144115188075855872, exceeding 4993/10000 by exactly
84177298230069/90071992547409920000. The output records every interval
bound and hashes both its input certificate and checker.

Unlike the original decimal-grid square-root calculation, the independent
checker uses 1<=v=1+U(a^2+2b^2)<4, starts with upper root 2, and performs
eight rational Newton steps,
rounding each upward to a 48-bit binary grid. Rational perfect squares are
handled exactly. Each final upper root is checked by squaring it. Each
integrand lower value is then rounded downward to the same grid. No
floating-point arithmetic or source-reported interval bound enters acceptance.
The Newton enclosure follows from
  ((u+v/u)/2)^2-v = (u^2-v)^2/(4u^2) >= 0.
Explicit exceptions preserve proof-relevant checks under Python optimization.

Tests cover the zero-tilt normalization, impossible strengthened targets,
geometric identities and other-sector exclusions at finite sample points,
invalid coefficient hypotheses, square-root enclosures, malformed rational
input, and missing/overlapping/degenerate parameter intervals. Finite point
controls supplement, and do not replace, the written continuum proof.

## Consequence and limits

The complete proof consists of Sections 1–4 together with the exact
finite checks. It shows that a normalized saving greater than 0.0007
below 1/2 is impossible within this fixed family, including arbitrary
positive shared even radii and all tilts in [0,1]. In particular, a
saving of 0.001 requires leaving at least one of these restrictions.
It neither establishes global optimality nor improves the explicit cover.
