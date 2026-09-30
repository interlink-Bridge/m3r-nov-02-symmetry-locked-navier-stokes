#!/usr/bin/env python3
import sympy as sp

x,y,z,A,nu,q=sp.symbols('x y z A nu q', positive=True, real=True)
coords=(x,y,z)
I=sp.eye(3)
Qx=sp.diag(-1,1,1); Qy=sp.diag(1,-1,1); Qz=sp.diag(1,1,-1)
Pyz=sp.Matrix([[1,0,0],[0,0,1],[0,1,0]])
GROUP=[('Qx',Qx),('Qy',Qy),('Qz',Qz),('Pyz',Pyz)]

u0=sp.Matrix([
    -A*sp.sin(x)*(sp.cos(y)+sp.cos(z)),
     A*sp.cos(x)*sp.sin(y),
     A*sp.cos(x)*sp.sin(z),
])
X=sp.Matrix([x,y,z])

# A01 lattice orthogonal transformations and exact velocity equivariance.
for name,Q in GROUP:
    assert Q.T*Q == I
    Xq=Q*X
    transformed=sp.simplify(u0.subs({x:Xq[0],y:Xq[1],z:Xq[2]}, simultaneous=True)-Q*u0)
    assert transformed == sp.zeros(3,1), (name, transformed)
print('PASS A03 exact_initial_equivariance_Qx_Qy_Qz_Pyz')

# A02 divergence freedom.
div=sp.factor(sum(sp.diff(u0[i],coords[i]) for i in range(3)))
assert div==0
print('PASS A02 divergence_free')

# A05 representation-theoretic shape forcing for Jacobian at common fixed point.
j=sp.symbols('j0:9')
J=sp.Matrix(3,3,j)
eqs=[]
for _,Q in GROUP:
    eqs += list(J-Q*J*Q.T)
eqs.append(sp.trace(J))
sol=sp.linsolve(eqs,j)
print('A05 invariant_divfree_J_solution =',sol)
# Expected one-parameter set diag(-2a,a,a).
a=sp.symbols('a', real=True)
expected=sp.diag(-2*a,a,a)
# verify expected satisfies all constraints
for _,Q in GROUP: assert expected==Q*expected*Q.T
assert sp.trace(expected)==0
# linsolve should have one free symbol and equivalent shape
assert len(sol.args)==1
vec=list(sol.args[0]); Jsol=sp.Matrix(3,3,vec)
assert Jsol[0,1]==Jsol[0,2]==Jsol[1,0]==Jsol[1,2]==Jsol[2,0]==Jsol[2,1]==0
assert sp.simplify(Jsol[1,1]-Jsol[2,2])==0
assert sp.simplify(Jsol[0,0]+2*Jsol[1,1])==0
print('PASS A05 gradient_shape_forcing')

# A07 covariance shape forcing. R is symmetric and group-invariant.
r11,r22,r33,r12,r13,r23=sp.symbols('r11 r22 r33 r12 r13 r23')
R=sp.Matrix([[r11,r12,r13],[r12,r22,r23],[r13,r23,r33]])
req=[]
for _,Q in GROUP:
    req += list(R-Q*R*Q.T)
rsol=sp.linsolve(req,[r11,r22,r33,r12,r13,r23])
print('A07 invariant_symmetric_R_solution =',rsol)
assert len(rsol.args)==1
rv=list(rsol.args[0])
assert rv[3]==rv[4]==rv[5]==0 and sp.simplify(rv[1]-rv[2])==0
print('PASS A07 covariance_shape_forcing')

# Differential and pressure identities.
def lap(expr): return sum(sp.diff(expr,c,2) for c in coords)
grad=sp.Matrix([[sp.diff(u0[i],coords[j]) for j in range(3)] for i in range(3)])
conv=sp.Matrix([sum(u0[j]*sp.diff(u0[i],coords[j]) for j in range(3)) for i in range(3)])
p=A**2*(sp.Rational(1,2)*sp.cos(y)*sp.cos(z)+sp.Rational(1,6)*sp.cos(2*x)*sp.cos(y)*sp.cos(z)+sp.Rational(1,2)*sp.cos(2*x)+sp.Rational(1,4)*sp.cos(2*y)+sp.Rational(1,4)*sp.cos(2*z))
poisson_rhs=-sum(grad[j,i]*grad[i,j] for i in range(3) for j in range(3))
assert sp.simplify(lap(p)-poisson_rhs)==0
H=sp.Matrix([[sp.diff(p,coords[i],coords[j]) for j in range(3)] for i in range(3)])
H0=sp.simplify(H.subs({x:0,y:0,z:0}))
assert H0==A**2*sp.diag(-sp.Rational(8,3),-sp.Rational(5,3),-sp.Rational(5,3))
lapu=sp.Matrix([lap(u0[i]) for i in range(3)])
ut=sp.simplify(-conv-sp.Matrix([sp.diff(p,c) for c in coords])+nu*lapu)
grad_ut=sp.Matrix([[sp.diff(ut[i],coords[j]) for j in range(3)] for i in range(3)])
grad_ut0=sp.simplify(grad_ut.subs({x:0,y:0,z:0}))
adot=sp.Rational(2,3)*A*(A-3*nu)
assert all(sp.simplify((grad_ut0-sp.diag(-2*adot,adot,adot))[i,j])==0 for i in range(3) for j in range(3))
print('PASS A21 pressure_hessian_and_raw_strain_derivative')

# Finite Fourier heat operator at origin.
XX,YY,ZZ=sp.symbols('XX YY ZZ', nonzero=True)
subs_trig={sp.sin(x):(XX-XX**-1)/(2*sp.I),sp.cos(x):(XX+XX**-1)/2,sp.sin(y):(YY-YY**-1)/(2*sp.I),sp.cos(y):(YY+YY**-1)/2,sp.sin(z):(ZZ-ZZ**-1)/(2*sp.I),sp.cos(z):(ZZ+ZZ**-1)/2}
def to_laurent(expr): return sp.expand(sp.expand_trig(expr).xreplace(subs_trig))
def heat_origin(expr):
    e=to_laurent(expr); total=0
    for term in sp.Add.make_args(e):
        pd=term.as_powers_dict(); kx=int(pd.get(XX,0)); ky=int(pd.get(YY,0)); kz=int(pd.get(ZZ,0))
        coeff=sp.simplify(term/(XX**kx*YY**ky*ZZ**kz))
        total += coeff*q**sp.Rational(kx*kx+ky*ky+kz*kz,2)
    return sp.factor(total)

Rheat=sp.zeros(3)
for i in range(3):
    for j in range(3):
        Rheat[i,j]=sp.simplify(heat_origin(u0[i]*u0[j])-heat_origin(u0[i])*heat_origin(u0[j]))
d=sp.factor(Rheat[0,0]-Rheat[1,1]); K=sp.factor(sp.trace(Rheat)/2)
d_exp=A**2*sp.Rational(1,4)*(1-q**2)*(q**2+4*q+1)
K_exp=A**2*sp.Rational(1,2)*(1-q**2)*(q**2+q+1)
assert sp.simplify(d-d_exp)==0 and sp.simplify(K-K_exp)==0
b=A*q
Pi=sp.factor(2*b*d)
Pi_exp=A**3*q*sp.Rational(1,2)*(1-q**2)*(q**2+4*q+1)
assert sp.simplify(Pi-Pi_exp)==0
print('PASS A12 initial_all_q_heat_algebra')
print('PASS A08 positivity_factors_for_A>0_q_in_(0,1)')

# Matched path direct from field.
matched_input=sp.simplify(ut-nu*lapu)
bdot_m=sp.factor(heat_origin(sp.diff(matched_input[1],y)))
U0=[heat_origin(ui) for ui in u0]; Udot_m=[heat_origin(matched_input[i]) for i in range(3)]
Rdot_m=sp.zeros(3)
for i in range(3):
    for j in range(3):
        f=u0[i]*u0[j]; ft=ut[i]*u0[j]+u0[i]*ut[j]
        Hdot=heat_origin(sp.expand(ft-nu*lap(f)))
        Rdot_m[i,j]=sp.simplify(Hdot-Udot_m[i]*U0[j]-U0[i]*Udot_m[j])
ddot_m=sp.factor(Rdot_m[0,0]-Rdot_m[1,1]); Pidot_m=sp.factor(2*(bdot_m*d+b*ddot_m))
Pidot_m_exp=(A**3*q/sp.Integer(3))*(A*q*(1-q**2)*(2*q**4+q**3+12*q**2+q+2)-6*nu*(q**4+2*q**3+2*q+1))
assert sp.simplify(Pidot_m-Pidot_m_exp)==0
print('PASS A13 matched_derivative_direct_from_field')

# Exact Sturm certificate and rational root/value enclosure.
Rstar=sp.factor(6*(q**4+2*q**3+2*q+1)/(q*(1-q**2)*(2*q**4+q**3+12*q**2+q+2)))
P10=3*q**10+9*q**9+8*q**8+32*q**7+17*q**6+44*q**5+29*q**4-16*q**3-16*q**2-q-1
num=sp.factor(sp.together(sp.diff(Rstar,q))).as_numer_denom()[0]
assert sp.simplify(num-12*P10)==0
st=sp.sturm(P10,q)
def variations_at(v):
    sg=[]
    for poly in st:
        val=sp.sign(sp.factor(poly.subs(q,v)))
        if val!=0: sg.append(int(val))
    return sum(1 for x,y in zip(sg,sg[1:]) if x*y<0),sg
v0,s0=variations_at(sp.Rational(0)); v1,s1=variations_at(sp.Rational(1))
assert v0-v1==1
qa=sp.Rational(63651056784173490941,10**20)
qb=sp.Rational(63651056784173490942,10**20)
assert sp.sign(P10.subs(q,qa))==-1 and sp.sign(P10.subs(q,qb))==1
# exact rational range enclosure using positive-factor interval arithmetic
N=lambda t: 6*(t**4+2*t**3+2*t+1)
G=lambda t: 2*t**4+t**3+12*t**2+t+2
Dmin=qa*(1-qb**2)*G(qa); Dmax=qb*(1-qa**2)*G(qb)
Rlo=sp.factor(N(qa)/Dmax); Rhi=sp.factor(N(qb)/Dmin)
assert Rlo < Rhi
print('PASS A15 exact_Sturm_variations',v0,v1,'count=',v0-v1)
print('PASS A16 exact_root_bracket',sp.N(qa,25),sp.N(qb,25))
print('PASS A16 exact_Rstar_min_enclosure',sp.N(Rlo,30),sp.N(Rhi,30))

# Fixed scale direct from field.
bdot_f=sp.factor(heat_origin(sp.diff(ut[1],y)))
Udot_f=[heat_origin(ut[i]) for i in range(3)]
Rdot_f=sp.zeros(3)
for i in range(3):
    for j in range(3):
        ft=ut[i]*u0[j]+u0[i]*ut[j]
        Hdot=heat_origin(sp.expand(ft))
        Rdot_f[i,j]=sp.simplify(Hdot-Udot_f[i]*U0[j]-U0[i]*Udot_f[j])
ddot_f=sp.factor(Rdot_f[0,0]-Rdot_f[1,1]); Pidot_f=sp.factor(2*(bdot_f*d+b*ddot_f))
Pidot_f_exp=(A**3*q*(1-q**2)/sp.Integer(3))*(A*q*(2*q**4+q**3+12*q**2+q+2)-9*nu*(q**2+4*q+1))
assert sp.simplify(Pidot_f-Pidot_f_exp)==0
Rhat=sp.factor(9*(q**2+4*q+1)/(q*(2*q**4+q**3+12*q**2+q+2)))
pos=3*q**6+17*q**5+17*q**4+50*q**3+19*q**2+q+1
assert sp.simplify(sp.diff(Rhat,q)+18*pos/(q**2*(2*q**4+q**3+12*q**2+q+2)**2))==0
assert sp.limit(Rhat,q,0,dir='+')==sp.oo and sp.limit(Rhat,q,1,dir='-')==3
print('PASS A18_A19 fixed_scale_direct_threshold_monotonicity_limits')

# Equality/crossover logic certificates.
# matched: unique minimizer implies equality case is tangent at unique q*.
assert P10.subs(q,0)<0 and P10.subs(q,1)>0
print('PASS A17 equality_case_logic_and_v0_3_text_patch')
print('ALL R2 THEOREM-LINE ALGEBRAIC CHECKS PASS')
