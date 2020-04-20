# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 11:39:00 2018

@author: Ilias
"""
import sys


def comInLoads(function, L, q, a=0, b=0):
    '''Combines aLL the functions to compute the reactions. Arg function is
    the appropriate function for the specific Load and is passed here through
    the GUI. Arg q can be either a vector or a tupLe(for the triagLe Load)'''
    if function[0]=='c':
        m1,m2,q1,q2 = getattr(sys.modules[__name__], function)(L,q, a, b)
        return m1,m2,q1,q2
    elif function[0]=='a':
        n1, n2 = getattr(sys.modules[__name__], function)(L,q, a, b)
        return n1, n2

def cFE(L, q,  a, b):
    Q1 = Q2 = q*L/2
    M1 = M2 = q*L*L/12
    return M1, -M2, -Q1, -Q2

def cPE(L, q,  a, b):
    Q1 = q*(L**4-(b**3)*(2*L-b)-a*(2*L**3-2*L*a**2+a**3))/(2*L**3)
    Q2 = q*(L**4-(a**3)*(2*L-a)-b*(2*L**3-2*L*b**2+b**3))/(2*L**3)
    M1 = q*(L**4-(b**3)*(4*L-3*b)-(a**2)*(6*L**2-8*a*L+3*a**2))/12*L**2
    M2 = q*(L**4-(a**3)*(4*L-3*a)-(b**2)*(6*L**2-8*b*L+3*b**2))/12*L**2
    return M1, -M2, -Q1, -Q2 


def cFT(L, q,  a, b):
    #gia ayjanomeno trigwniko. gia fortio pou meiwnetai anapoda ta Q1,Q2,M1,M2
    q1, q2 = q
    Q1 = 3*q2*L/20
    Q2 = 7*q2*L/20
    M1 = q2*L*L/30
    M2 = q2*L*L/20

    if q1 == 0:
        return M1, -M2, -Q1, -Q2

    else:
        Q1, Q2 = Q2, Q1
        M1, M2 = M2, M1
        return M1, -M2, -Q1, -Q2



def cPT(L, q,  a, b):
    #gia ayjanomeno trigwniko. gia fortio pou meiwnetai anapoda ta Q1,Q2,M1,M2
    b=L-b
    q1, q2 = q
    Q1 = -a**2*q2/(2*L) + b**2*q2/(2*L) + 3*a**4*q2/(4*L**3) - 3*b**4*q2/(4*L**3) - 2*a**5*q2/(5*L**4) + 2*b**5*q2/(5*L**4)
    Q2 = -3*a**4*q2/(4*L**3) + 3*b**4*q2/(4*L**3) + 2*a**5*q2/(5*L**4) - 2*b**5*q2/(5*L**4)
    M1 = -a**3*q2/(3*L) + b**3*q2/(3*L) + a**4*q2/(2*L**2) - b**4*q2/(2*L**2) - a**5*q2/(5*L**3) + b**5*q2/(5*L**3)
    M2 = -a**4*q2/(4*L**2) + b**4*q2/(4*L**2) + a**5*q2/(5*L**3) - b**5*q2/(5*L**3)
    if q1 == 0:
        return M1, -M2, -Q1, -Q2 
    else: 
        Q1, Q2 = Q2, Q1
        M1, M2 = M2, M1
        return M1, -M2, -Q1, -Q2 



def cPo(L, q,  a, b):
    Q1 = (q*(b**2)/(L**2))*(3-2*b/L)
    Q2 = (q*(a**2)/(L**2))*(3-2*a/L)
    M1 = q*a*(b**2)/(L**2)
    M2 = q*b*(a**2)/(L**2)
    return M1, -M2, -Q1, -Q2




def cM(L, q,  a, b):
    Q1 = 6*q*a*b/(L**3)
    Q2 = -6*q*a*b/(L**3)
    M1 = q*b*(2-3*b/L)/L
    M2 = q*a*(2-3*a/L)/L
    return M1, M2, -Q1, Q2



def aFE(L, q,  a, b):
    N1 = q*(L-a-b)*(L-a+b)/(2*L)
    N2 = q*(L-a-b)*(L+a-b)/(2*L)
    return -N1, -N2



def aPE(L, q,  a, b):
    N1 = q*(L-a-b)*(L-a+b)/(2*L)
    N2 = q*(L-a-b)*(L+a-b)/(2*L)
    return -N1, -N2




def aFT(L, q,  a, b):
    if q[0]==0:q=q[1]
    else:q=q[0]
    N1 = q*(L-a-b)*(L-a+b)/(2*L)
    N2 = q*(L-a-b)*(L+a-b)/(2*L)
    return -N1/2, -N2/2




def aPT(L, q,  a, b):
    if q[0]==0:q=q[1]
    else:q=q[0]
    N1 = q*(L-a-b)*(L-a+b)/(2*L)
    N2 = q*(L-a-b)*(L+a-b)/(2*L)
    return -N1/2, -N2/2



def aPo(L, q,  a, b):
    N1 = q*b/L
    N2 = q*a/L
    return -N1, -N2




def aM(L, q,  a, b): 
    M1 = q*b/L
    M2 = q*a/L
    return -M1, -M2

