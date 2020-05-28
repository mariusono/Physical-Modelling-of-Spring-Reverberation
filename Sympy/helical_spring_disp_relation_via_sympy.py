# -*- coding: utf-8 -*-
"""
Created on Fri May 29 00:03:04 2020

@author: MSOI
"""



from sympy import *
import numpy as np
import matplotlib.pyplot as plt
import re
import copy


I,eta,Kappa,k,Dxx,u,Dxxxx,q_slash,q,gamma,alpha,sig_t,sig_l,Dxmin,Dxpl,zeta,theta,Mx,sig_0,sig_1 = symbols('I,eta,Kappa,k,Dxx,u,Dxxxx,q_slash,q,gamma,alpha,sig_t,sig_l,Dxmin,Dxpl,zeta,theta,Mx,sig_0,sig_1')

U,Z,beta = symbols('U,Z,beta')
z,w,x,t,beta,K,k,mu,h = symbols('z,w,x,t,beta,K,k,mu,h')
from sympy import I


u_lm0_nm0,u_lm0_nm1,u_lm0_nm2,u_lm1_nm0,u_lm1_nm1,u_lm1_nm2,u_lm2_nm0,u_lm2_nm1,u_lm2_nm2,u_lp0_np0,u_lp0_np1,u_lp0_np2,u_lp1_np0,u_lp1_np1,u_lp1_np2,u_lp2_np0,u_lp2_np1,u_lp2_np2,zeta_lm0_nm0,zeta_lm0_nm1,zeta_lm0_nm2,zeta_lm1_nm0,zeta_lm1_nm1,zeta_lm1_nm2,zeta_lm2_nm0,zeta_lm2_nm1,zeta_lm2_nm2,zeta_lp0_np0,zeta_lp0_np1,zeta_lp0_np2,zeta_lp1_np0,zeta_lp1_np1,zeta_lp1_np2,zeta_lp2_np0,zeta_lp2_np1,zeta_lp2_np2 = symbols('u_lm0_nm0,u_lm0_nm1,u_lm0_nm2,u_lm1_nm0,u_lm1_nm1,u_lm1_nm2,u_lm2_nm0,u_lm2_nm1,u_lm2_nm2,u_lp0_np0,u_lp0_np1,u_lp0_np2,u_lp1_np0,u_lp1_np1,u_lp1_np2,u_lp2_np0,u_lp2_np1,u_lp2_np2,zeta_lm0_nm0,zeta_lm0_nm1,zeta_lm0_nm2,zeta_lm1_nm0,zeta_lm1_nm1,zeta_lm1_nm2,zeta_lm2_nm0,zeta_lm2_nm1,zeta_lm2_nm2,zeta_lp0_np0,zeta_lp0_np1,zeta_lp0_np2,zeta_lp1_np0,zeta_lp1_np1,zeta_lp1_np2,zeta_lp2_np0,zeta_lp2_np1,zeta_lp2_np2')

u_lm0_np0,u_lm0_np1,u_lm0_np2,u_lm1_np0,u_lm1_np1,u_lm1_np2,u_lm2_np0,u_lm2_np1,u_lm2_np2,u_lp0_nm0,u_lp0_nm1,u_lp0_nm2,u_lp1_nm0,u_lp1_nm1,u_lp1_nm2,u_lp2_nm0,u_lp2_nm1,u_lp2_nm2,zeta_lm0_np0,zeta_lm0_np1,zeta_lm0_np2,zeta_lm1_np0,zeta_lm1_np1,zeta_lm1_np2,zeta_lm2_np0,zeta_lm2_np1,zeta_lm2_np2,zeta_lp0_nm0,zeta_lp0_nm1,zeta_lp0_nm2,zeta_lp1_nm0,zeta_lp1_nm1,zeta_lp1_nm2,zeta_lp2_nm0,zeta_lp2_nm1,zeta_lp2_nm2 = symbols('u_lm0_np0,u_lm0_np1,u_lm0_np2,u_lm1_np0,u_lm1_np1,u_lm1_np2,u_lm2_np0,u_lm2_np1,u_lm2_np2,u_lp0_nm0,u_lp0_nm1,u_lp0_nm2,u_lp1_nm0,u_lp1_nm1,u_lp1_nm2,u_lp2_nm0,u_lp2_nm1,u_lp2_nm2,zeta_lm0_np0,zeta_lm0_np1,zeta_lm0_np2,zeta_lm1_np0,zeta_lm1_np1,zeta_lm1_np2,zeta_lm2_np0,zeta_lm2_np1,zeta_lm2_np2,zeta_lp0_nm0,zeta_lp0_nm1,zeta_lp0_nm2,zeta_lp1_nm0,zeta_lp1_nm1,zeta_lp1_nm2,zeta_lp2_nm0,zeta_lp2_nm1,zeta_lp2_nm2')


## Dispersion relation helical spring

U,Z,beta = symbols('U,Z,beta')
from sympy import I

expr1 = diff(U*exp(I*(w*t+beta*x)),t,t) + K**2*( diff(U*exp(I*(w*t+beta*x)),x,x,x,x) + 2*q**2*diff(U*exp(I*(w*t+beta*x)),x,x) + q**4*U*exp(I*(w*t+beta*x)) ) - q**2*gamma**2*(diff(Z*exp(I*(w*t+beta*x)),x)-U*exp(I*(w*t+beta*x)))
expr2 = diff(Z*exp(I*(w*t+beta*x)),t,t) - gamma**2*(diff(Z*exp(I*(w*t+beta*x)),x,x)-diff(U*exp(I*(w*t+beta*x)),x))

expr = expr1

exprNew = expr.expand(U)
factors_U_expr1 = collect(exprNew, U, evaluate=False)

for key in list(factors_U_expr1.keys()):
    factors_U_expr1[key] = simplify(factors_U_expr1[key])


expr = factors_U_expr1[1]

exprNew = expr.expand(Z)
factors_Z_expr1 = collect(exprNew, Z, evaluate=False)

for key in list(factors_Z_expr1.keys()):
    factors_Z_expr1[key] = simplify(factors_Z_expr1[key])


expr = expr2

exprNew = expr.expand(U)
factors_U_expr2 = collect(exprNew, U, evaluate=False)

for key in list(factors_U_expr2.keys()):
    factors_U_expr2[key] = simplify(factors_U_expr2[key])


expr = factors_U_expr2[1]

exprNew = expr.expand(Z)
factors_Z_expr2 = collect(exprNew, Z, evaluate=False)

for key in list(factors_Z_expr2.keys()):
    factors_Z_expr2[key] = simplify(factors_Z_expr2[key])




ana = Matrix([[factors_U_expr1[U],factors_Z_expr1[Z]],
              [factors_U_expr2[U],factors_Z_expr2[Z]]]);


solsNew = solve(ana.det(),w)



# Plot the model dispersion
beta_num = np.arange(0,2000+10,10)

k_val = 1/44100
K_val = 0.06
q_val = 835
gamma_val = 1980
eta_val = 0.411
theta_val = 0.000411
alpha_val = 0.5


w_num_1 = np.zeros(np.shape(beta_num))
w_num_3 = np.zeros(np.shape(beta_num))
for iBeta in range(0,len(beta_num)):
    w_num_1[iBeta] = solsNew[1].evalf(subs={K:K_val,gamma:gamma_val,q:q_val,beta:beta_num[iBeta]})
    w_num_3[iBeta] = solsNew[3].evalf(subs={K:K_val,gamma:gamma_val,q:q_val,beta:beta_num[iBeta]})


plt.figure()
plt.plot(beta_num,w_num_1/(2*np.pi),'-',linewidth=2)
plt.plot(beta_num,w_num_3/(2*np.pi),'-',linewidth=2)
plt.grid()
plt.show()


plt.figure()
plt.plot(beta_num,w_num_1/(2*np.pi),'-',linewidth=2)
plt.grid()
plt.ylim([0,20000])
plt.show()




## Numerical Disp relation helical spring implicit fd scheme: 
## no loss terms in the following expr:

from sympy import I

l,n,z = symbols('l,n,z')

expr_u = ( (-2*u_lm0_nm0 + u_lm0_nm1 + u_lm0_np1)/k**2 + eta*K*k*(4*u_lm0_nm0 - 2*u_lm0_nm1 - 2*u_lm0_np1 - 2*u_lm1_nm0 + u_lm1_nm1 + u_lm1_np1 - 2*u_lp1_nm0 + u_lp1_nm1 + u_lp1_np1)/(k**2*h**2) ) \
        + K**2*( (6*u_lm0_nm0 - 4*u_lm1_nm0 + u_lm2_nm0 - 4*u_lp1_nm0 + u_lp2_nm0)/h**4 + 2*q_slash**2*((-2*u_lm0_nm0 + u_lm1_nm0 + u_lp1_nm0)/h**2) + q_slash**4*u_lm0_nm0 ) \
        - ( gamma**2*q**2*alpha*((zeta_lm0_nm0 - zeta_lm1_nm0)/h) - gamma**2*q**2*alpha*u_lm0_nm0 + gamma**2*q**2*(1-alpha)*(1/2)*((zeta_lm0_np1 - zeta_lm1_np1)/h + (zeta_lm0_nm1 - zeta_lm1_nm1)/h) + gamma**2*q**2*(1-alpha)*(1/2)*(u_lm0_np1 + u_lm0_nm1) )

expr_u = expr_u.subs(u_lm0_nm0,U*z**0*exp(0*I*beta*h))
expr_u = expr_u.subs(u_lm0_nm1,U*z**-1*exp(0*I*beta*h))
expr_u = expr_u.subs(u_lm0_np1,U*z**1*exp(0*I*beta*h))

expr_u = expr_u.subs(u_lm1_nm0,U*z**0*exp(-1*I*beta*h))
expr_u = expr_u.subs(u_lm1_nm1,U*z**-1*exp(-1*I*beta*h))
expr_u = expr_u.subs(u_lm1_np1,U*z**1*exp(-1*I*beta*h))

expr_u = expr_u.subs(u_lp1_nm0,U*z**0*exp(1*I*beta*h))
expr_u = expr_u.subs(u_lp1_nm1,U*z**-1*exp(1*I*beta*h))
expr_u = expr_u.subs(u_lp1_np1,U*z**1*exp(1*I*beta*h))

expr_u = expr_u.subs(u_lm2_nm0,U*z**0*exp(-2*I*beta*h))
expr_u = expr_u.subs(u_lm2_nm1,U*z**-1*exp(-2*I*beta*h))
expr_u = expr_u.subs(u_lm2_np1,U*z**1*exp(-2*I*beta*h))

expr_u = expr_u.subs(u_lp2_nm0,U*z**0*exp(2*I*beta*h))
expr_u = expr_u.subs(u_lp2_nm1,U*z**-1*exp(2*I*beta*h))
expr_u = expr_u.subs(u_lp2_np1,U*z**1*exp(2*I*beta*h))

expr_u = expr_u.subs(zeta_lm0_nm0,Z*z**0*exp(0*I*beta*h))
expr_u = expr_u.subs(zeta_lm0_nm1,Z*z**-1*exp(0*I*beta*h))
expr_u = expr_u.subs(zeta_lm0_np1,Z*z**1*exp(0*I*beta*h))

expr_u = expr_u.subs(zeta_lm1_nm0,Z*z**0*exp(-1*I*beta*h))
expr_u = expr_u.subs(zeta_lm1_nm1,Z*z**-1*exp(-1*I*beta*h))
expr_u = expr_u.subs(zeta_lm1_np1,Z*z**1*exp(-1*I*beta*h))

expr_u = expr_u.subs(zeta_lp1_nm0,Z*z**0*exp(1*I*beta*h))
expr_u = expr_u.subs(zeta_lp1_nm1,Z*z**-1*exp(1*I*beta*h))
expr_u = expr_u.subs(zeta_lp1_np1,Z*z**1*exp(1*I*beta*h))


expr = expr_u

exprNew = expr.expand(U)
factors_U_expr1 = collect(exprNew, U, evaluate=False)

# Latex printing
init_printing(pretty_print=True)

# Disable Latex printing
init_printing(pretty_print=False)

for key in list(factors_U_expr1.keys()):
    factors_U_expr1[key] = simplify(factors_U_expr1[key])


expr = factors_U_expr1[1]

exprNew = expr.expand(Z)
factors_Z_expr1 = collect(exprNew, Z, evaluate=False)

for key in list(factors_Z_expr1.keys()):
    factors_Z_expr1[key] = simplify(factors_Z_expr1[key])


expr_zeta = ( (-2*zeta_lm0_nm0 + zeta_lm0_nm1 + zeta_lm0_np1)/k**2 + theta*gamma**2*k**2*(4*zeta_lm0_nm0 - 2*zeta_lm0_nm1 - 2*zeta_lm0_np1 - 2*zeta_lm1_nm0 + zeta_lm1_nm1 + zeta_lm1_np1 - 2*zeta_lp1_nm0 + zeta_lp1_nm1 + zeta_lp1_np1)/(k**2*h**2) ) \
        - ( gamma**2*alpha*((-2*zeta_lm0_nm0 + zeta_lm1_nm0 + zeta_lp1_nm0)/h**2) - gamma**2*alpha*(-u_lm0_nm0 + u_lp1_nm0)/h + gamma**2*(1-alpha)*(1/2)*((-2*zeta_lm0_np1 + zeta_lm1_np1 + zeta_lp1_np1)/h**2 + (-2*zeta_lm0_nm1 + zeta_lm1_nm1 + zeta_lp1_nm1)/h**2) - gamma**2*(1-alpha)*(1/2)*((-u_lm0_np1 + u_lp1_np1)/h + (-u_lm0_nm1 + u_lp1_nm1)/h) )



expr_zeta = expr_zeta.subs(zeta_lm0_nm0,Z*z**0*exp(0*I*beta*h))
expr_zeta = expr_zeta.subs(zeta_lm0_nm1,Z*z**-1*exp(0*I*beta*h))
expr_zeta = expr_zeta.subs(zeta_lm0_np1,Z*z**1*exp(0*I*beta*h))

expr_zeta = expr_zeta.subs(zeta_lm1_nm0,Z*z**0*exp(-1*I*beta*h))
expr_zeta = expr_zeta.subs(zeta_lm1_nm1,Z*z**-1*exp(-1*I*beta*h))
expr_zeta = expr_zeta.subs(zeta_lm1_np1,Z*z**1*exp(-1*I*beta*h))

expr_zeta = expr_zeta.subs(zeta_lp1_nm0,Z*z**0*exp(I*beta*h))
expr_zeta = expr_zeta.subs(zeta_lp1_nm1,Z*z**-1*exp(I*beta*h))
expr_zeta = expr_zeta.subs(zeta_lp1_np1,Z*z**1*exp(I*beta*h))

expr_zeta = expr_zeta.subs(u_lm0_nm0,U*z**0*exp(0*I*beta*h))
expr_zeta = expr_zeta.subs(u_lm0_nm1,U*z**-1*exp(0*I*beta*h))
expr_zeta = expr_zeta.subs(u_lm0_np1,U*z**1*exp(0*I*beta*h))

expr_zeta = expr_zeta.subs(u_lm1_nm0,U*z**0*exp(-1*I*beta*h))
expr_zeta = expr_zeta.subs(u_lm1_nm1,U*z**-1*exp(-1*I*beta*h))
expr_zeta = expr_zeta.subs(u_lm1_np1,U*z**1*exp(-1*I*beta*h))

expr_zeta = expr_zeta.subs(u_lp1_nm0,U*z**0*exp(1*I*beta*h))
expr_zeta = expr_zeta.subs(u_lp1_nm1,U*z**-1*exp(1*I*beta*h))
expr_zeta = expr_zeta.subs(u_lp1_np1,U*z**1*exp(1*I*beta*h))


expr = expr_zeta

exprNew = expr.expand(U)
factors_U_expr2 = collect(exprNew, U, evaluate=False)

for key in list(factors_U_expr2.keys()):
    factors_U_expr2[key] = simplify(factors_U_expr2[key])


expr = factors_U_expr2[1]

exprNew = expr.expand(Z)
factors_Z_expr2 = collect(exprNew, Z, evaluate=False)

for key in list(factors_Z_expr2.keys()):
    factors_Z_expr2[key] = simplify(factors_Z_expr2[key])


ana = Matrix([[factors_U_expr1[U],factors_Z_expr1[Z]],
              [factors_U_expr2[U],factors_Z_expr2[Z]]]);


expr = ana.det()
expr = expr.subs(z,exp(I*w*k)) # expression in Matlab script main_generate_spring_irs_and_everything_else.m line 169

# Try to solve it numerically.. done in Matlab

