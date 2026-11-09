#!/usr/bin/env python3
"""Exact cycle-q2 forcing reduction and discriminant-residue checks.

Provenance: extends the connection reduction by benediamond in
https://github.com/hodge-computer/sextics/pull/223#issuecomment-6712406587
and fqlx's unfinished cycle current in PR #221. Input is the protected
cycle-q2 instance, on X2=1. This is partial progress, not a certificate.
Requires Python 3 and SymPy; writes JSON beside this file.
"""
import json
from pathlib import Path
import sympy as s

# Independent reconstruction of the universal degree-six residue connection.
y,c,f,r=s.symbols('y c f r')
K=s.QQ.frac_field(c,f)
F=y**6+c*y+f
g=s.diff(F,y)
Delta=6**6*f**5-5**5*c**6
Fin=s.Poly(F,y,domain=K)
def red(v):
    return s.Poly(v,y,domain=K).rem(Fin).as_expr()
def canon(v):
    p=s.Poly(red(v),y,domain=K)
    return s.Add(*(s.factor(a)*y**i for (i,),a in p.terms()))
H=s.invert(g,F,y)
assert red(g*H)==1

def delta(v,t):
    return canon(s.diff(v,t)-(y if t==c else 1)*H*s.diff(v,y))
basis=[y,y**2,y**3,y**4,y**5+5*c/s.Integer(6)]
def coords(v):
    p=s.Poly(red(v),y,domain=K)
    u=s.Matrix([s.factor(p.nth(i)) for i in range(1,6)])
    assert s.cancel(p.nth(0)-5*c*u[4]/6)==0
    return u
Bc=s.Matrix.hstack(*(coords(delta(e,c)) for e in basis))
Bf=s.Matrix.hstack(*(coords(delta(e,f)) for e in basis))
assert s.cancel(s.trace(Bc)-s.diff(Delta,c)/(2*Delta))==0
assert s.cancel(s.trace(Bf)-s.diff(Delta,f)/(2*Delta))==0
assert all(s.cancel(v)==0 for v in Bf.diff(c)-Bc.diff(f)+Bc*Bf-Bf*Bc)
assert red(delta(y,c)+y*H)==0

# At a generic discriminant point the critical root r is nonzero.
critical={c:-6*r**5,f:5*r**6}
assert s.expand(Delta.subs(critical))==0
v=s.Matrix([r**(-i) for i in range(5)])
w=s.Matrix([s.Rational(i+1,30)*r**i for i in range(5)])
R=v*w.T
for t,B in [(c,Bc),(f,Bf)]:
    residue=(B*Delta).applyfunc(s.cancel).subs(critical)/s.diff(Delta,t).subs(critical)
    assert all(s.cancel(a)==0 for a in residue-R)
assert R.rank()==1 and s.trace(R)==s.Rational(1,2)
assert all(s.cancel(a)==0 for a in R*R-R/2)
phi=coords(y*H)
phi_residue=(phi*Delta).applyfunc(s.cancel).subs(critical)/s.diff(Delta,f).subs(critical)
assert all(s.cancel(a)==0 for a in phi_residue+r*v/30)
assert all(s.cancel(a)==0 for a in (s.eye(5)-2*R)*phi_residue)

# Physical source binding. All four physical base coordinates stay independent.
x,p,z,q=s.symbols('x p z q')
base=(x,p,z,q) # (X0,X1,Y0,Y1), with X2=1 and y=Y2.
c0=(x*z+p**3*q)/2
f0=x**6+p**6+1+z**6+q**6+x**3*p*z*q/2
a=x**2*p**2*z**2*q**2
F0=y**6+c0*y+f0
g0=s.diff(F0,y)
# Bind to the protected record when this packet is run inside its repository.
# The record has one-line exact F and A entries; no YAML dependency is needed.
record=Path(__file__).resolve().parents[3]/'data/cycle-q2/instance.rec'
if record.is_file():
    entries=dict(line.split(': ',1) for line in record.read_text().splitlines()
                 if line.startswith(('F: ','A: ','q: ')))
    original=s.symbols('X0 X1 X2 Y0 Y1 Y2')
    chart=dict(zip(original,(x,p,s.Integer(1),z,q,y)))
    assert s.expand(s.sympify(entries['F']).subs(chart)-F0)==0
    assert s.expand(s.sympify(entries['A']).subs(chart)-a*y)==0
    assert entries['q']=='2'
T=lambda u:s.expand(x*s.diff(u,x)-z*s.diff(u,z))
h=T(f0)
h1=T(h)
h2=T(h1)
lam=s.expand(3*h1**2-h*h2)
assert T(c0)==T(a)==0
assert s.diff(h,y)==0 and h!=0
fx,fz=s.diff(f0,x),s.diff(f0,z)
# This is the useful Jacobian identity for the same conical denominator h.
assert s.expand(s.diff(c0,x)*fz-s.diff(c0,z)*fx+h/2)==0
b=s.Matrix([-2*fz/h,0,2*fx/h,0])
assert s.cancel(sum(bj*s.diff(c0,t) for bj,t in zip(b,base))-1)==0
assert s.cancel(sum(bj*s.diff(f0,t) for bj,t in zip(b,base)))==0
sigma=a*lam/(2*h**4)
# kappa=div(sigma*b) has one explicit numerator over h^6.
Al=s.expand(a*lam)
P=s.expand(h*(-s.diff(Al,x)*fz+s.diff(Al,z)*fx)
           +5*Al*(fz*s.diff(h,x)-fx*s.diff(h,z)))
kappa=P/h**6
J=sigma*b*y
assert s.cancel(sum(s.diff(sigma*bj,t) for bj,t in zip(b,base))-kappa)==0
# A literal rational identity, even before reduction modulo F0:
# sum_j delta_j J_j = kappa*y - sigma*y/g0.
# Coefficient derivatives, including all denominator derivatives, are retained.
transgression=sum(s.diff(Jj,t)-s.diff(F0,t)*s.diff(Jj,y)/g0
                  for Jj,t in zip(J,base))
assert s.cancel(transgression-kappa*y+sigma*y/g0)==0
assert s.Poly(P,base).length()==15
assert s.Poly(lam,base).length()==6

out={
 'status':'partial progress only; no rational solution of the remaining PDE supplied',
 'provenance':'https://github.com/hodge-computer/sextics/pull/223#issuecomment-6712406587',
 'chart':'X2=1; x=X0, p=X1, z=Y0, q=Y1, y=Y2',
 'c':str(c0),'f':str(f0),'a':str(a),'h':str(h),'lambda':str(lam),
 'sigma':str(sigma),'b_c':[str(s.factor(vv)) for vv in b],
 'new_forcing_numerator_P':str(s.factor(P)),
 'new_forcing_denominator':'h**6',
 'new_forcing':'-(P/h**6) * (1,0,0,0,0)^T',
 'current_translation':'u_j -> u_j - sigma*b_c[j]*(1,0,0,0,0)^T',
 'identity':'sum_j delta_j(sigma*b_c[j]*y) = (P/h**6)*y - sigma*y/g',
 'P_monomials':15,
 'nonzero_test':{'point':{'x':1,'p':1,'z':1,'q':1},
                 'h':str(h.subs(dict.fromkeys(base,1))),
                 'P':str(P.subs(dict.fromkeys(base,1)))},
 'discriminant':'Delta=46656*f**5-3125*c**6',
 'critical_parameterization':'c=-6*r**5, f=5*r**6, r!=0',
 'residue_outer_vectors':{'v':[str(vv) for vv in v], 'w':[str(vv) for vv in w]},
 'residue_identities':['R=v*w^T','rank R=1','trace R=1/2','R^2=R/2',
                       'trace B_t = (partial_t Delta)/(2 Delta)',
                       'residue_f(y/g) = -r*v/30'],
 'checks':['universal connection reconstructed independently','connection flatness',
           'both transverse residues','forcing residue projector',
           'exact physical coefficient binding','physical b_c(c)=1 and b_c(f)=0',
           'literal rational transgression identity','all coefficient and denominator derivatives']
}
if record.is_file():
    out['checks'].append('exact protected instance.rec F, A and pole order binding')
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print('PASS connection flatness, logarithmic trace, rank-one half-residue, physical Jacobian,')
print('     and exact one-component forcing reduction with all four base variables retained.')
print('P has',out['P_monomials'],'monomials; h(1,1,1,1)=',out['nonzero_test']['h'],
      '; P(1,1,1,1)=',out['nonzero_test']['P'])
