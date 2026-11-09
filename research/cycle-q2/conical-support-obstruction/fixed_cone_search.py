#!/usr/bin/env python3
"""Exact complete fixed-D=h^m homogeneous search for cycle-q2.

If a dual certificate is exported, it is checked directly over ZZ; no rank
inference from finite characteristic is used. Fixed-D completeness follows
from monic division by F^2 in Y2 and diagonal character averaging, as explained
in the accompanying report. This does not exclude other denominator supports.
"""
import argparse,json,time,hashlib,ast
from fractions import Fraction
from pathlib import Path
from itertools import product
start=time.time()
P=argparse.ArgumentParser();P.add_argument('--m',type=int,default=0);P.add_argument('--verify',type=Path);P.add_argument('--instance',type=Path,default=Path(__file__).with_name('instance.rec'));args=P.parse_args();m=args.m
assert m>=0
zero=(0,)*6
couplings=[(3,1,0,1,1,0),(0,3,1,0,1,1),(1,0,3,1,0,1)]
A=(2,2,3,2,2,1)
F2={tuple(6*(i==j) for j in range(6)):2 for i in range(6)}
F2.update({v:1 for v in couplings})
h={(6,0,0,0,0,0):6,(0,0,0,6,0,0):-6,(3,1,0,1,1,0):1}

def addexp(a,b):return tuple(x+y for x,y in zip(a,b))
def mul(a,b):
 c={}
 for e,v in a.items():
  for f,w in b.items():
   t=addexp(e,f);c[t]=c.get(t,0)+v*w
 return {e:v for e,v in c.items() if v}
def power(p,n):
 r={zero:1}
 for _ in range(n):r=mul(r,p)
 return r
def deriv(p,i):
 out={}
 for e,v in p.items():
  if e[i]:
   t=list(e);t[i]-=1;out[tuple(t)]=v*e[i]
 return out
def comps(n,k):
 if k==1:
  yield(n,);return
 for i in range(n+1):
  for r in comps(n-i,k-1):yield(i,)+r
# Exact binding to the supplied protected instance snapshot (stdlib only).
varnames=['X0','X1','X2','Y0','Y1','Y2']
def parsepoly(node):
 if isinstance(node,ast.Expression):return parsepoly(node.body)
 if isinstance(node,ast.Constant) and type(node.value) is int:return {zero:Fraction(node.value)}
 if isinstance(node,ast.Name) and node.id in varnames:
  return {tuple(int(j==varnames.index(node.id)) for j in range(6)):Fraction(1)}
 if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):return {e:-v for e,v in parsepoly(node.operand).items()}
 if isinstance(node,ast.BinOp):
  if isinstance(node.op,ast.Pow):
   assert isinstance(node.right,ast.Constant) and type(node.right.value) is int and node.right.value>=0
   return power(parsepoly(node.left),node.right.value)
  l,r=parsepoly(node.left),parsepoly(node.right)
  if isinstance(node.op,ast.Mult):return mul(l,r)
  if isinstance(node.op,ast.Div):
   assert set(r)=={zero} and r[zero]
   return {e:v/r[zero] for e,v in l.items()}
  if isinstance(node.op,(ast.Add,ast.Sub)):
   out=dict(l)
   for e,v in r.items():out[e]=out.get(e,0)+(v if isinstance(node.op,ast.Add) else -v)
   return {e:v for e,v in out.items() if v}
 raise ValueError('Unsupported expression in instance record')
record=args.instance.read_bytes()
fields=dict(line.split(':',1) for line in record.decode().splitlines() if line and not line.startswith('#') and ':' in line)
assert fields['id'].strip()=='cycle-q2' and fields['q'].strip()=='2'
assert {e:2*v for e,v in parsepoly(ast.parse(fields['F'].strip(),mode='eval')).items()}==F2
assert parsepoly(ast.parse(fields['A'].strip(),mode='eval'))=={A:1}
rowspace={tuple(sum(n[i]*couplings[i][j] for i in range(3))%6 for j in range(6)) for n in product(range(6),repeat=3)}
G=[w for w in product(range(6),repeat=6) if all(sum(a*b for a,b in zip(w,v))%6==0 for v in couplings)]
assert len(G)==432 and len(rowspace)==108
assert len(G)*len(rowspace)==6**6
assert all(sum(a*b for a,b in zip(g,v))%6==0 for g in G for v in rowspace)
unknowns=[]
for e in comps(6*m+7,6):
 if e[5]>=12:continue
 for i in range(6):
  if tuple((e[j]-A[j]-(i==j))%6 for j in range(6)) in rowspace:unknowns.append((i,e))
print('m',m,'unknowns',len(unknowns),'elapsed',round(time.time()-start,2),flush=True)
if m:
 pref=[{e:2*v for e,v in mul(h,deriv(F2,i)).items()} for i in range(6)]
 divpref=mul(F2,h)
 hpref=[{e:m*v for e,v in mul(F2,deriv(h,i)).items()} for i in range(6)]
 rhs={addexp(e,A):2*v for e,v in power(h,m+1).items()}
else:
 pref=[{e:2*v for e,v in deriv(F2,i).items()} for i in range(6)]
 divpref=F2;hpref=[{}]*6;rhs={A:2}
columns=[];rowset=set(rhs)
for i,e in unknowns:
 col={}
 for p in [pref[i],hpref[i]]:
  for f,v in p.items():
   t=addexp(e,f);col[t]=col.get(t,0)+v
 if e[i]:
  em=list(e);em[i]-=1
  for f,v in divpref.items():
   t=addexp(em,f);col[t]=col.get(t,0)-e[i]*v
 col={t:v for t,v in col.items() if v}
 rowset.update(col);columns.append(col)
rows=sorted(rowset);index={e:i for i,e in enumerate(rows)}
print('rows',len(rows),'nnz',sum(map(len,columns)),'elapsed',round(time.time()-start,2),flush=True)
if args.verify:
 saved=json.loads(args.verify.read_text())
 assert saved['m']==m and saved['unknowns']==len(unknowns) and saved['rows']==len(rows)
 relation=saved['exact_dual_certificate']
 cert={tuple(e):int(v) for e,v in relation['functional']}
 assert len(cert)==len(relation['functional']) and set(cert)<=rowset
 assert all(sum(cert.get(e,0)*v for e,v in col.items())==0 for col in columns)
 pairing=sum(cert.get(e,0)*v for e,v in rhs.items())
 assert pairing==int(relation['rhs_pairing']) and pairing!=0
 print('PASS independent exact ZZ dual verification; annihilates all',len(columns),
       'full-space columns and pairs nontrivially with source; elapsed',round(time.time()-start,2))
 raise SystemExit(0)
from flint import fmpz_mat
# Integer transpose: each row is one unknown numerator monomial's image.
AT=fmpz_mat(len(columns),len(rows))
for j,col in enumerate(columns):
 for e,v in col.items():AT[j,index[e]]=v
print('integer nullspace starts',flush=True)
Z,nullity=AT.nullspace()
print('nullity',nullity,'elapsed',round(time.time()-start,2),flush=True)
found=None
for k in range(nullity):
 dot=sum(int(Z[index[e],k])*v for e,v in rhs.items())
 if dot:
  cert={e:int(Z[j,k]) for j,e in enumerate(rows) if Z[j,k]}
  # Independently verify every coefficient of the exact dual relation.
  assert all(sum(cert.get(e,0)*v for e,v in col.items())==0 for col in columns)
  assert sum(cert.get(e,0)*v for e,v in rhs.items())==dot
  found={'functional':[[list(e),str(v)] for e,v in cert.items()], 'rhs_pairing':str(dot)}
  break
out={'instance':'cycle-q2','m':m,'D':'h**m','normal_degree_bound':'deg_Y2 N_i < 12',
     'numerator_degree':6*m+7,'group_order':len(G),'group_rowspace_order':len(rowspace),
     'unknowns':len(unknowns),'rows':len(rows),'nnz':sum(map(len,columns)),
     'nullity_transpose':nullity,'exact_dual_certificate':found,
     'arithmetic':'integer matrices using python-flint; exact integer dual checked against every column',
     'result':'inconsistent over QQ (and CC)' if found else 'no obstruction found in computed left nullspace',
     'elapsed_seconds':time.time()-start}
path=Path(__file__).with_name(f'fixed_cone_m{m}.json')
path.write_text(json.dumps(out,indent=2)+'\n')
print(out['result'],'dual terms',len(found['functional']) if found else 0,'elapsed',round(time.time()-start,2),flush=True)
