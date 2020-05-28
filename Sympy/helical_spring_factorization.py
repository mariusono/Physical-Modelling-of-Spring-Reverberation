# -*- coding: utf-8 -*-
"""
Created on Thu May 28 23:54:40 2020

@author: MSOI
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 23:47:59 2020

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

expr1 = (I+eta*Kappa*k*Dxx)*(1/k**2)*(u**3-2*u**2+u) + Kappa**2*(Dxxxx*u**2+2*q_slash**2*Dxx*u**2+I*q_slash**4*u**2) - ( q**2*gamma**2*alpha*Dxmin*zeta**2 - q**2*gamma**2*alpha*u**2 + q**2*gamma**2*(1-alpha)*Dxmin*(1/(2))*(zeta**3+zeta) - q**2*gamma**2*(1-alpha)*(1/(2))*(u**3+u) ) + sig_t/k*(u**3-u) 
expr2 =  (I+theta*gamma**2*k**2*Dxx)*(1/k**2)*(zeta**3-2*zeta**2+zeta) - ( + gamma**2*alpha*Dxx*zeta**2 - gamma**2*alpha*Dxpl*u**2 + gamma**2*(1-alpha)*Dxx*(1/(2))*(zeta**3+zeta) - gamma**2*(1-alpha)*Dxpl*(1/(2))*(u**3+u) - (sig_l/k)*(zeta**3-zeta) )               
   


expr = expr1

exprNew = expr.expand(u)
factors = collect(exprNew, u, evaluate=False)

for key in list(factors.keys()):
    factors[key] = simplify(factors[key])

factors

expr = factors[1]

exprNew = expr.expand(zeta)
factors_zeta = collect(exprNew, zeta, evaluate=False)

for key in list(factors_zeta.keys()):
    factors_zeta[key] = simplify(factors_zeta[key])

factors_zeta




