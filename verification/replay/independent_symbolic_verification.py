#!/usr/bin/env python3
"""Independent symbolic replay for M3R-NOV-02 v0.3 referee-patched candidate.

Reconstructs the explicit initial field and checks both the calorically matched
path and the fixed-heat-scale adversarial control directly from the field.
SymPy is the only non-stdlib dependency.
"""
import sympy as sp

x,y,z,A,nu,q=sp.symbols('x y z A nu q', positive=True, real=True)
coords=(x,y,z)
u=sp.Matrix([
    -A*sp.sin(x)*(sp.cos(y)+sp.cos(z)),
    A*sp.cos(x)*sp.sin(y),
    A*sp.cos(x)*sp.sin(z),
])

def lap(expr):
    return sum(sp.diff(expr,c,2) for c in coords)

# Core differential identities
div=sp.simplify(sum(sp.diff(u[i],coords[i]) for i in range(3)))
assert div == 0

grad=sp.Matrix([[sp.diff(u[i],coords[j]) for j in range(3)] for i in range(3)])
conv=sp.Matrix([sum(u[j]*sp.diff(u[i],coords[j]) for j in range(3)) for i in range(3)])

p=A**2*(
    sp.Rational(1,2)*sp.cos(y)*sp.cos(z)
    +sp.Rational(1,6)*sp.cos(2*x)*sp.cos(y)*sp.cos(z)
    +sp.Rational(1,2)*sp.cos(2*x)
    +sp.Rational(1,4)*sp.cos(2*y)
    +sp.Rational(1,4)*sp.cos(2*z)
)
poisson_rhs=-sum(grad[j,i]*grad[i,j] for i in range(3) for j in range(3))
assert sp.simplify(lap(p)-poisson_rhs) == 0

H=sp.Matrix([[sp.diff(p,coords[i],coords[j]) for j in range(3)] for i in range(3)])
H0=sp.simplify(H.subs({x:0,y:0,z:0}))
H0_expected=A**2*sp.diag(-sp.Rational(8,3),-sp.Rational(5,3),-sp.Rational(5,3))
assert H0 == H0_expected

lapu=sp.Matrix([lap(u[i]) for i in range(3)])
ut=sp.simplify(-conv-sp.Matrix([sp.diff(p,c) for c in coords])+nu*lapu)
grad_ut=sp.Matrix([[sp.diff(ut[i],coords[j]) for j in range(3)] for i in range(3)])
grad_ut0=sp.simplify(grad_ut.subs({x:0,y:0,z:0}))
adot=sp.Rational(2,3)*A*(A-3*nu)
expected_grad_ut0 = sp.diag(-2*adot,adot,adot)
assert all(sp.simplify(grad_ut0[i,j]-expected_grad_ut0[i,j]) == 0 for i in range(3) for j in range(3))

# Exact heat semigroup at the origin via finite Laurent/Fourier algebra.
X,Y,Z=sp.symbols('X Y Z', nonzero=True)
subs_trig={
    sp.sin(x):(X-X**-1)/(2*sp.I), sp.cos(x):(X+X**-1)/2,
    sp.sin(y):(Y-Y**-1)/(2*sp.I), sp.cos(y):(Y+Y**-1)/2,
    sp.sin(z):(Z-Z**-1)/(2*sp.I), sp.cos(z):(Z+Z**-1)/2,
}

def to_laurent(expr):
    return sp.expand(sp.expand_trig(expr).xreplace(subs_trig))

def heat_origin(expr):
    e=to_laurent(expr)
    total=0
    for term in sp.Add.make_args(e):
        pd=term.as_powers_dict()
        kx=int(pd.get(X,0)); ky=int(pd.get(Y,0)); kz=int(pd.get(Z,0))
        coeff=sp.simplify(term/(X**kx*Y**ky*Z**kz))
        total += coeff*q**sp.Rational(kx*kx+ky*ky+kz*kz,2)
    return sp.factor(total)

R=sp.zeros(3)
for i in range(3):
    for j in range(3):
        R[i,j]=sp.simplify(heat_origin(u[i]*u[j])-heat_origin(u[i])*heat_origin(u[j]))

d=sp.factor(R[0,0]-R[1,1])
K=sp.factor(sp.trace(R)/2)
d_expected=A**2*sp.Rational(1,4)*(1-q**2)*(q**2+4*q+1)
K_expected=A**2*sp.Rational(1,2)*(1-q**2)*(q**2+q+1)
assert sp.simplify(d-d_expected) == 0
assert sp.simplify(K-K_expected) == 0

# Matched physical-time derivative with tau'(t)=-1.
matched_input=sp.simplify(ut-nu*lapu)
bdot=sp.factor(heat_origin(sp.diff(matched_input[1],y)))
bdot_expected=sp.Rational(2,3)*A**2*q**3
assert sp.simplify(bdot-bdot_expected) == 0

U0=[heat_origin(ui) for ui in u]
Udot0=[heat_origin(matched_input[i]) for i in range(3)]
Rdot=sp.zeros(3)
for i in range(3):
    for j in range(3):
        f=u[i]*u[j]
        ft=ut[i]*u[j]+u[i]*ut[j]
        Hdot=heat_origin(sp.expand(ft-nu*lap(f)))
        Rdot[i,j]=sp.simplify(Hdot-Udot0[i]*U0[j]-U0[i]*Udot0[j])

ddot=sp.factor(Rdot[0,0]-Rdot[1,1])
ddot_expected=(A**3/sp.Integer(3))*q*(1-q**2)*(q**4+4*q**2+1)-nu*A**2*(q**4+2*q**3+2*q+1)
assert sp.simplify(ddot-ddot_expected) == 0

b=A*q
Pidot=sp.factor(2*(bdot*d+b*ddot))
Pidot_expected=(A**3*q/sp.Integer(3))*(A*q*(1-q**2)*(2*q**4+q**3+12*q**2+q+2)-6*nu*(q**4+2*q**3+2*q+1))
assert sp.simplify(Pidot-Pidot_expected) == 0

Rstar=sp.factor(6*(q**4+2*q**3+2*q+1)/(q*(1-q**2)*(2*q**4+q**3+12*q**2+q+2)))
num=sp.factor(sp.together(sp.diff(Rstar,q))).as_numer_denom()[0]
P10=3*q**10+9*q**9+8*q**8+32*q**7+17*q**6+44*q**5+29*q**4-16*q**3-16*q**2-q-1
assert sp.simplify(num-12*P10) == 0
root_count=sp.polys.polytools.count_roots(P10,0,1)
assert root_count == 1
roots=[r for r in sp.nroots(P10,n=30,maxsteps=200) if abs(complex(r).imag)<1e-20 and 0<float(sp.re(r))<1]
assert len(roots)==1
qstar=sp.N(sp.re(roots[0]),20)
rmin=sp.N(Rstar.subs(q,qstar),20)

# Fixed-heat-scale physical-time derivative: tau is held constant.
bdot_fixed=sp.factor(heat_origin(sp.diff(ut[1],y)))
bdot_fixed_expected=sp.Rational(2,3)*A*q*(A*q**2-3*nu)
assert sp.simplify(bdot_fixed-bdot_fixed_expected) == 0

Udot_fixed=[heat_origin(ut[i]) for i in range(3)]
Rdot_fixed=sp.zeros(3)
for i in range(3):
    for j in range(3):
        ft=ut[i]*u[j]+u[i]*ut[j]
        Hdot=heat_origin(sp.expand(ft))
        Rdot_fixed[i,j]=sp.simplify(Hdot-Udot_fixed[i]*U0[j]-U0[i]*Udot_fixed[j])

ddot_fixed=sp.factor(Rdot_fixed[0,0]-Rdot_fixed[1,1])
ddot_fixed_expected=(A**2*(1-q**2)/sp.Integer(3))*(A*q*(q**4+4*q**2+1)-3*nu*(q**2+4*q+1))
assert sp.simplify(ddot_fixed-ddot_fixed_expected) == 0

Pidot_fixed=sp.factor(2*(bdot_fixed*d+b*ddot_fixed))
Pidot_fixed_expected=(A**3*q*(1-q**2)/sp.Integer(3))*(A*q*(2*q**4+q**3+12*q**2+q+2)-9*nu*(q**2+4*q+1))
assert sp.simplify(Pidot_fixed-Pidot_fixed_expected) == 0

Rhat=sp.factor(9*(q**2+4*q+1)/(q*(2*q**4+q**3+12*q**2+q+2)))
Rhat_der=sp.factor(sp.diff(Rhat,q))
positive_poly=3*q**6+17*q**5+17*q**4+50*q**3+19*q**2+q+1
Rhat_der_expected=-18*positive_poly/(q**2*(2*q**4+q**3+12*q**2+q+2)**2)
assert sp.simplify(Rhat_der-Rhat_der_expected) == 0
assert sp.limit(Rhat,q,0,dir='+') == sp.oo
assert sp.limit(Rhat,q,1,dir='-') == 3

print('PASS divergence_free')
print('PASS pressure_poisson')
print('PASS pressure_hessian_origin', H0)
print('PASS raw_strain_derivative adot =', sp.factor(adot))
print('PASS heat_covariance d =', d)
print('PASS heat_covariance K =', K)
print('PASS matched bdot =', bdot)
print('PASS matched ddot =', ddot)
print('PASS matched Pidot =', Pidot)
print('PASS matched_threshold_derivative_polynomial =', P10)
print('PASS matched_root_count_(0,1) =', root_count)
print('q_star =', qstar)
print('R_star_min =', rmin)
print('PASS fixed_scale bdot =', bdot_fixed)
print('PASS fixed_scale ddot =', ddot_fixed)
print('PASS fixed_scale Pidot =', Pidot_fixed)
print('PASS fixed_scale Rhat =', Rhat)
print('PASS fixed_scale Rhat_derivative =', Rhat_der)
print('PASS fixed_scale limits = +infinity, 3')
print('ALL MATCHED AND FIXED-SCALE SYMBOLIC CHECKS PASS')
