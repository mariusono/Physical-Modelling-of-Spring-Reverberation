# -*- coding: utf-8 -*-
"""
Created on Fri May 29 00:07:48 2020

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

## SLINKY MODEL

## Model dispersion:

sigma_0, sigma_1 = symbols('sigma_0,sigma_1')

expr = diff(U*exp(I*(w*t+beta*x)),t,t) + K**2*( diff(U*exp(I*(w*t+beta*x)),x,x,x,x) ) + 2*sigma_0*diff(U*exp(I*(w*t+beta*x)),t) - 2*sigma_1*diff(U*exp(I*(w*t+beta*x)),t,x,x)

expr = expr.subs(U,1)

#expr_toSolve = expr.rewrite(sin)
#expr.rewrite(sin)

sols = solve(expr,w)



