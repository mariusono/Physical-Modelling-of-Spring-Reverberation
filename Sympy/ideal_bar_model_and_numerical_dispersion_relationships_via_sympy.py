# -*- coding: utf-8 -*-
"""
Created on Fri May 29 00:06:29 2020

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

## Test on 1d wave

expr = diff(U*exp(I*(w*t+beta*x)),t,t) - gamma**2*( diff(U*exp(I*(w*t+beta*x)),x,x) )

expr = expr.subs(U,1)
#expr.rewrite(sin)

sols = solve(expr,w)

## Test on the ideal bar

## Analytical:

expr = diff(U*exp(I*(w*t+beta*x)),t,t) + K**2*( diff(U*exp(I*(w*t+beta*x)),x,x,x,x) )

expr = expr.subs(U,1)
#expr.rewrite(sin)

sols = solve(expr,w)



expr = diff(U*exp(I*(w*t+beta*x)),t,t) + K**2*( diff(U*exp(I*(w*t+beta*x)),x,x,x,x) )
expr = expr.subs(U,1)
#expr.rewrite(sin)

sols = solve(expr,w)


## Explitict FDS: 
expr_u = u_lm0_np1 - ( (2-6*mu**2)*u_lm0_nm0 + 4*mu**2*(u_lp1_nm0+u_lm1_nm0) - mu**2*(u_lp2_nm0+u_lm2_nm0) - u_lm0_nm1 )

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


expr_u = expr_u.subs(U,1)

expr_u.rewrite(sin)

expr_toSolve = expr_u.rewrite(sin)
expr_toSolve = expr_toSolve.subs(z,exp(I*w*k))

#expr_toSolve.evalf(subs={K:K_val,mu:mu_val,k:k_val,h:h_val,beta:50,w:380*2*np.pi})

expr_toSolve = expr_toSolve.rewrite(sin)

sols = solve(expr_toSolve,w)



## Implcit FDS:

expr_u = theta*(1/k**2)*(u_lm0_np1 - 2*u_lm0_nm0 + u_lm0_nm1) \
        + (1-theta)*0.5*((1/k**2)*(u_lp1_np1-2*u_lp1_nm0+u_lp1_nm1)+(1/k**2)*(u_lm1_np1-2*u_lm1_nm0+u_lm1_nm1)) \
        + K**2*(1/h**4)*(u_lp2_nm0 - 4*u_lp1_nm0 + 6*u_lm0_nm0 - 4*u_lm1_nm0 + u_lm2_nm0)
     


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

expr_u = expr_u.subs(U,1)

#expr_u = expr_u.subs(z,exp(I*w*k))
#sols = solve(expr_u,w) # not working ? 

# # Better:
expr_u.rewrite(sin)

expr_toSolve = expr_u.rewrite(sin)
#expr_toSolve = expr_u
expr_toSolve = expr_toSolve.subs(z,exp(I*w*k))

#expr_toSolve.evalf(subs={K:K_val,mu:mu_val,k:k_val,h:h_val,beta:50,w:380*2*np.pi})

expr_toSolve = expr_toSolve.rewrite(sin)

sols = solve(expr_toSolve,w)


