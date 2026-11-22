# Local cube lines at the four pure endpoints

An invariant line has coordinates $(a,b,c,d)=(X,rX,Y,sY)$, with the alternate standard chart at $r=\infty$ or $s=\infty$. Write $f(r)=1+r^6+\kappa r+2\kappa r^5$, $g(s)=1+s^6+\kappa s+3\kappa s^5$, and $\kappa=5/7$. If one ratio is $0$ or $\infty$, a perfect cube requires the complementary ratio to be a root of its binary sextic, giving 24 geometric points. All the roots are simple and nonzero. The local calculation uses the full four-dimensional Grassmann chart.

At an $X$-endpoint use $c=Y$, $d=wX+(s+z)Y$ and $g(s)=0$. At zero use $(a,b)=(X,uX+vY)$; at infinity use $(a,b)=(uX+vY,X)$. At a $Y$-endpoint use $a=Y$, $b=wX+(r+z)Y$ and $f(r)=0$, with $(c,d)$ given by those respective two choices. In all four cases the central restricted form is $X^6$. Taking $T=Y/X$, normalize the central value of its quadratic cube root $Q$ to $1$. The three coefficients of $Q$ are locally determined by rows $T^0,T^1,T^2$ of $H-Q^3$, with diagonal linear coefficient $-3$.

At $r=0$, the remaining rows $T^3,T^4,T^5,T^6$ have linear parts $\kappa s v,2\kappa s u,2\kappa s v+g'(s)w,g'(s)z$. At $s=0$ they have linear parts $2\kappa r v,\kappa r u,\kappa r v+f'(r)w,f'(r)z$. Thus the 12 zero endpoints are reduced isolated points.

At infinity, the three rows $T^4,T^5,T^6$ have linear parts $C u,C v+D w,D z$, where $(C,D)=(2\kappa s,g'(s))$ at $r=\infty$ and $(\kappa r,f'(r))$ at $s=\infty$. The six indicated rows solve for $u,w,z$ and all three coefficients of $Q$ uniquely as formal functions of $v$. The $T^3$ row has zero linear term; it is odd in $v$ by sign equivariance, since $v,w$ change sign and $T$ is negated. The remaining test is its cubic coefficient.

The exact checker works over the whole complementary sextic algebra modulo $11$ and computes through $v^3$. The pivots are units in that algebra. With $s$ the complementary ratio, the obstruction coefficients are

| Infinite endpoint | Complementary binary mod 11 | Coefficient of $v^3$ in the residual $T^3$ row |
|---|---|---|
| $X$ | $s^6-s^5-4s+1$ | $-3s^5-s^4-5s^2+s$ |
| $Y$ | $s^6+3s^5-4s+1$ | $-3s^5-4s^4-4s^3+3s^2+4s$ |

Each coefficient has gcd one with its respective binary. The same formal elimination is therefore integral at 11, and the characteristic-zero norm of each coefficient is nonzero. Each of the 12 infinite endpoints has local length three and is isolated. The check concerns these 24 endpoints only. The other invariant and non-invariant lines require the separate projective calculation. [Arithmetic checker](pure_lines.py)
