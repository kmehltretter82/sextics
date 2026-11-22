"""Check the twelve lifts of the nine-parameter quotient on the smooth mod-11 sample."""
from itertools import product
from fiber_current import AE,BO,QBASE,AB,table,ev,deriv,add
from fiber_point import xbase

p=11
class E:
    """a+b*z, z^2+z+1=0, over F_11."""
    def __init__(self,a=0,b=0):self.a,self.b=a%p,b%p
    def __add__(self,o):
        if not isinstance(o,E):o=E(o)
        return E(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return E(-self.a,-self.b)
    def __sub__(self,o):return self+-o
    def __rsub__(self,o):return -self+o
    def __mul__(self,o):
        if not isinstance(o,E):o=E(o)
        return E(self.a*o.a-self.b*o.b,self.a*o.b+self.b*o.a-self.b*o.b)
    __rmul__=__mul__
    def __pow__(self,n):
        if n<0:return self**(n%120)
        if not n:return E(1)
        return (self*self)**(n//2)*(self if n%2 else 1)
    def __truediv__(self,o):
        if not isinstance(o,E):o=E(o)
        assert o.a or o.b
        return self*o**119
    def __eq__(self,o):
        if not isinstance(o,E):o=E(o)
        return self.a==o.a and self.b==o.b

def main():
    q=table(QBASE+[AB],[4,4,5,10,7,1]);Q=table(QBASE,xbase[5:]);A=table(AE,[5,3,1,8,9,5,6,4]);B=table(BO,[7,9,5,2,4,5,0,6])
    qt=table(QBASE,[4,6,7,3,9]);At=table(AE,[0,3,5,0,8,6,2,0]);Bt=table(BO,[2,5,5,5,3,0,8,0])
    r,s=4,2;Z=3;L=8;z=E(0,1)
    assert (z**2+z+1)==0 and z!=1
    Qhat={m:int(v)*pow(L,-1,p)%p for m,v in Q.items()}
    for al in [6,5]:
        pt=[1,r,al,s*al]
        assert al*al%p==Z and ev(q,pt)%p==0 and ev(Qhat,pt)%p==1
        assert ev(A,pt)*ev(A,pt)%p==ev(B,pt)*ev(B,pt)%p
    H=(1+r**6+7*(r+2*r**5)+7*r*s*Z+14*r*s*Z**2+(1+s**6+7*(s+3*s**5))*Z**3)%p
    assert H==L**3%p==6
    lifts=set()
    for e,eps,sgn in product(range(3),[1,-1],[1,-1]):
        lam=z**e
        h=add((A,eps),(B,eps*sgn))
        dh=add((At,eps),(Bt,eps*sgn))
        alpha=(sgn*6)%p
        pt=[1,r,alpha,s*alpha]
        assert ev(q,pt)%p==ev(h,pt)%p==0
        conjugate=[1,r,-alpha,-s*alpha]
        assert ev(h,conjugate)%p
        mappedL=lam*ev(Q,pt)
        assert mappedL**3==H and all(lam*v/mappedL==Qhat[m] for m,v in Q.items())
        hvec=tuple(h.get(m,0) for m in AE+BO)
        lifts.add((e,hvec,alpha))
        def partial(f):
            b,c,d=[ev(deriv(f,i),pt) for i in [1,2,3]]
            return b%p,alpha*d%p,(c+s*d)%p
        qr,qs,qa=partial(q);hr,hs,ha=partial(h)
        br,bs,ba=[-lam*y for y in partial(Q)];beta=-mappedL
        jac=(ev(qt,pt)*ha-ev(dh,pt)*qa)%p
        er,es,ez=(qa*hr-ha*qr)%p,(qa*hs-ha*qs)%p,(qr*hs-qs*hr)%p
        ey=er*bs-es*br+ba*ez
        # Changing h by an overall sign changes both Dh and its residue
        # Jacobian by that sign. The variable involution sends R to -R;
        # scaling Q multiplies 3Q^2 and hence the residue Jacobian by lam².
        ell=E(sgn*8)/lam**2
        fy=1+s**6+7*(s+3*s**5);fys=6*s**5+7*(1+15*s**4)
        entries=[]
        for i,(v,pr,ps) in enumerate([(7*r*s,7*s,7*r),(14*r*s,14*s,14*r),(fy,0,fys)],1):
            entries.append(alpha*(3*v*beta*ey-beta**2*(er*ps-es*pr))-2*i*v*beta**2*ez)
        a,b,c=entries;b+=alpha**2*c;a+=alpha**2*b
        assert [4*ell*y/(3*jac) for y in (a,b,c)]==[9,9,6]
    assert len(lifts)==12
    print("12 distinct lifts over F_121, Qhat =",[Qhat[m] for m in QBASE],", Z =",Z,", L^3 =",H,", source for every lift (9,9,6)")
if __name__=="__main__":main()
