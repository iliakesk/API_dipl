#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 17:28:04 2020

@author: ilias
"""
import numpy as np
import ireact3d as ir3


class Load():
    def __init__(self, load_value, element_type, load_type, a, b):
        self.load_value = load_value
        self.load_type = load_type
        self.start = a
        self.end = b
        

#epistrefei ta epikomvia fortia poy prokaloyntai apo to antistoixo epi toy meloys
def getEdgeValues(axes, L, q, load_type, a=0, b=0):
    axx, axy, axz = axes
    if load_type=='M':
        Rx, Ry, Rz = getMomentReactions(L, axx, axy, axz, q, a, b, load_type)
    elif load_type[1]=='T':#triangular load
        Rx, Ry, Rz = getTriangReactions(L, axx, axy, axz, q, a, b, load_type)
    elif load_type[1]=='R':#trapezoidal load
        Rx, Ry, Rz = getTrapezoidReactions(L, axx, axy, axz, q, a, b, load_type)
    elif load_type[1]=='E':#uniform load
        Rx, Ry, Rz = getUniformReactions(L, axx, axy, axz, q, a, b, load_type)
    elif load_type=='Po':#uniform load
        b = L-a
        Rx, Ry, Rz = getPointReactions(L, axx, axy, axz, q, a, b, load_type)
    Nstart, Nend = Rx[0]*axx, Rx[1]*axx
    zMstart, zMend, yQstart, yQend = Ry[0]*axz, Ry[1]*axz, Ry[2]*axy, Ry[3]*axy
    yMstart, yMend, zQstart, zQend = Rz[0]*axy, Rz[1]*axy, Rz[2]*axz, Rz[3]*axz
    if load_type=='M':
        Ps = 0 + yQstart + zQstart
        Pe = 0 + yQend + zQend
        Ms = Nstart + yMstart + zMstart
        Me = Nend + yMend + zMend
    else:
        Ps = Nstart + yQstart + zQstart
        Pe = Nend + yQend + zQend
        Ms = 0 + yMstart + zMstart
        Me = 0 + yMend + zMend

    loadsstart = np.array([Ps[0], Ps[1], Ps[2], Ms[0], Ms[1], Ms[2]])
    loadsstop = np.array([Pe[0], Pe[1], Pe[2], Me[0], Me[1], Me[2]])
    return loadsstart, loadsstop


def getComponents(axx, axy, axz, q, m=False):
    #xrisimopoieitai mono gia to analysis, an to fortio einai simeiako i omoiomorfa katanemimeno
#    ikanopoieitai i prwti sinthiki giati einai tis morfis [Plx1,Ply1,Plz1,Plx2,Ply2,Plz2]
#    an einai trigwnoko ikanopoiitai i deyteri shnthhkh giati einai [0,[Plx1,Ply1,Plz1,Plx2,Ply2,Plz2]] h
#    [[Plx1,Ply1,Plz1,Plx2,Ply2,Plz2], 0]. i teleytaia synthkh an einai trapezoides
#    PREPEI NA SPASEI SE KATI PIO KATANOHTO
    if len(q)==3:
        P=q
        if np.dot(P,axx)==0:Plx=0
        else:Plx = np.dot(P,axx)*np.linalg.norm((np.dot(P,axx))*axx)/abs(np.dot(P,axx))
        if np.dot(P,axy)==0:Ply=0
        else:Ply = np.dot(P,axy)*np.linalg.norm((np.dot(P,axy))*axy)/abs(np.dot(P,axy))
        if np.dot(P,axz)==0:Plz=0
        else:Plz = np.dot(P,axz)*np.linalg.norm((np.dot(P,axz))*axz)/abs(np.dot(P,axz))
        Plx = 0 if abs(Plx)<0.5 else Plx
        Ply = 0 if abs(Ply)<0.5 else Ply
        Plz = 0 if abs(Plz)<0.5 else Plz
        if m:
            return 0,0,0,Plx,Ply,Plz
        else:
            return Plx,Ply,Plz,0,0,0
    elif q[0]==0:
        P=q[1]
        if np.dot(P,axx)==0:Plx=0
        else:Plx = np.dot(P,axx)*np.linalg.norm((np.dot(P,axx))*axx)/abs(np.dot(P,axx))
        if np.dot(P,axy)==0:Ply=0
        else:Ply = np.dot(P,axy)*np.linalg.norm((np.dot(P,axy))*axy)/abs(np.dot(P,axy))
        if np.dot(P,axz)==0:Plz=0
        else:Plz = np.dot(P,axz)*np.linalg.norm((np.dot(P,axz))*axz)/abs(np.dot(P,axz))
        Plx = 0 if abs(Plx)<0.5 else Plx
        Ply = 0 if abs(Ply)<0.5 else Ply
        Plz = 0 if abs(Plz)<0.5 else Plz
        return 0,0,0,Plx,Ply,Plz
    elif q[1]==0:
        P=q[0]
        if np.dot(P,axx)==0:Plx=0
        else:Plx = np.dot(P,axx)*np.linalg.norm((np.dot(P,axx))*axx)/abs(np.dot(P,axx))
        if np.dot(P,axy)==0:Ply=0
        else:Ply = np.dot(P,axy)*np.linalg.norm((np.dot(P,axy))*axy)/abs(np.dot(P,axy))
        if np.dot(P,axz)==0:Plz=0
        else:Plz = np.dot(P,axz)*np.linalg.norm((np.dot(P,axz))*axz)/abs(np.dot(P,axz))
        Plx = 0 if abs(Plx)<0.5 else Plx
        Ply = 0 if abs(Ply)<0.5 else Ply
        Plz = 0 if abs(Plz)<0.5 else Plz
        return Plx,Ply,Plz,0,0,0
    else:
        Plx1,Ply1,Plz1 = getComponents(axx, axy, axz, q[0])
        Plx2,Ply2,Plz2 = getComponents(axx, axy, axz, q[1])
        return Plx1,Ply1,Plz1,Plx2,Ply2,Plz2


def getMomentReactions(L, axx, axy, axz, q, a, b, load_type):
    P=q
    if np.dot(P,axx)==0:Plx=0
    else:Plx = np.dot(P,axx)*np.linalg.norm((np.dot(P,axx))*axx)/abs(np.dot(P,axx))
    if np.dot(P,axy)==0:Ply=0
    else:Ply = np.dot(P,axy)*np.linalg.norm((np.dot(P,axy))*axy)/abs(np.dot(P,axy))
    if np.dot(P,axz)==0:Plz=0
    else:Plz = np.dot(P,axz)*np.linalg.norm((np.dot(P,axz))*axz)/abs(np.dot(P,axz))
    Rx = ir3.comInLoads('a'+load_type, L, Plx, a, b)
    Ry = ir3.comInLoads('c'+load_type, L, Ply, a, b)
    Rz = ir3.comInLoads('c'+load_type, L, Plz, a, b)
    return Rx, Ry, Rz


def getPointReactions(L, axx, axy, axz, q, a, b, load_type):
    P=q
    if np.dot(P,axx)==0:Plx=0
    else:Plx = np.dot(P,axx)*np.linalg.norm((np.dot(P,axx))*axx)/abs(np.dot(P,axx))
    if np.dot(P,axy)==0:Ply=0
    else:Ply = np.dot(P,axy)*np.linalg.norm((np.dot(P,axy))*axy)/abs(np.dot(P,axy))
    if np.dot(P,axz)==0:Plz=0
    else:Plz = np.dot(P,axz)*np.linalg.norm((np.dot(P,axz))*axz)/abs(np.dot(P,axz))
    Rx = ir3.comInLoads('a'+load_type, L, Plx, a, b)
    Ry = ir3.comInLoads('c'+load_type, L, Ply, a, b)
    Rz = ir3.comInLoads('c'+load_type, L, Plz, a, b)
    return Rx, Ry, Rz

def getUniformReactions(L, axx, axy, axz, q, a, b, load_type):
    P=q
    if np.dot(P,axx)==0:Plx=0
    else:Plx = np.dot(P,axx)*np.linalg.norm((np.dot(P,axx))*axx)/abs(np.dot(P,axx))
    if np.dot(P,axy)==0:Ply=0
    else:Ply = np.dot(P,axy)*np.linalg.norm((np.dot(P,axy))*axy)/abs(np.dot(P,axy))
    if np.dot(P,axz)==0:Plz=0
    else:Plz = np.dot(P,axz)*np.linalg.norm((np.dot(P,axz))*axz)/abs(np.dot(P,axz))
    Rx = ir3.comInLoads('a'+load_type, L, Plx, a, b)
    Ry = ir3.comInLoads('c'+load_type, L, Ply, a, b)
    Rz = ir3.comInLoads('c'+load_type, L, Plz, a, b)
    return Rx, Ry, Rz

def getTriangReactions(L, axx, axy, axz, q, a, b, load_type):
    if q[0]==0:
        P=q[1]
        if np.dot(P,axx)==0:Plx=0
        else:Plx = np.dot(P,axx)*np.linalg.norm((np.dot(P,axx))*axx)/abs(np.dot(P,axx))
        if np.dot(P,axy)==0:Ply=0
        else:Ply = np.dot(P,axy)*np.linalg.norm((np.dot(P,axy))*axy)/abs(np.dot(P,axy))
        if np.dot(P,axz)==0:Plz=0
        else:Plz = np.dot(P,axz)*np.linalg.norm((np.dot(P,axz))*axz)/abs(np.dot(P,axz))
        Rx = ir3.comInLoads('a'+load_type, L, (0,Plx), a, b)
        Ry = ir3.comInLoads('c'+load_type, L, (0,Ply), a, b)
        Rz = ir3.comInLoads('c'+load_type, L, (0,Plz), a, b)     
    else:
        P=q[0]
        if np.dot(P,axx)==0:Plx=0
        else:Plx = np.dot(P,axx)*np.linalg.norm((np.dot(P,axx))*axx)/abs(np.dot(P,axx))
        if np.dot(P,axy)==0:Ply=0
        else:Ply = np.dot(P,axy)*np.linalg.norm((np.dot(P,axy))*axy)/abs(np.dot(P,axy))
        if np.dot(P,axz)==0:Plz=0
        else:Plz = np.dot(P,axz)*np.linalg.norm((np.dot(P,axz))*axz)/abs(np.dot(P,axz))
        Rx = ir3.comInLoads('a'+load_type, L, (Plx,0), a, b)
        Ry = ir3.comInLoads('c'+load_type, L, (Ply,0), a, b)
        Rz = ir3.comInLoads('c'+load_type, L, (Plz,0), a, b) 
    return Rx, Ry, Rz
        
        
def getTrapezoidReactions(L, axx, axy, axz, q, a, b, load_type):
    if (np.array(np.absolute(q[0]))<=np.array(np.absolute(q[1]))).all():
        q_t = np.array(q[1]) - np.array(q[0])
        q_t = q_t.tolist()
        q_t = (0,q_t)
        q_u = q[0]
        load_typeT = load_type[0]+'T'
        load_typeE = load_type[0]+'E'
        Rxt, Ryt, Rzt = getTriangReactions(L, axx, axy, axz, q_t, a, b, load_typeT)
        Rxe, Rye, Rze = getUniformReactions(L, axx, axy, axz, q_u, a, b, load_typeE)
        Rx = Rxt+Rxe
        Ry = Ryt+Rye
        Rz = Rzt+Rze
    elif (np.array(np.absolute(q[0]))>=np.array(np.absolute(q[1]))).all():
        q_t = np.array(q[0]) - np.array(q[1])
        q_t = q_t.tolist()
        q_t = (q_t,0)
        q_u = q[1]
        load_typeT = load_type[0]+'T'
        load_typeE = load_type[0]+'E'
        Rxt, Ryt, Rzt = getTriangReactions(L, axx, axy, axz, q_t, a, b, load_typeT)
        Rxe, Rye, Rze = getUniformReactions(L, axx, axy, axz, q_u, a, b, load_typeE)
        Rx = Rxt+Rxe
        Ry = Ryt+Rye
        Rz = Rzt+Rze
    return Rx, Ry, Rz



#
#nodal_loads = {}
#
#
#def beamLoad(axes, L, snode, enode, q, load_type, a=0, b=0):
#    '''Adds to the Loads dataframe the loads that are computed
#    by the inner reactions of the members'''
#    for node in [snode, enode]:
#        if node not in nodal_loads.keys():
#            nodal_loads[node]=...............np,array()
#        else:
#            nodal_loads[node] += ...............
#    loadsstart, loadsstop = getLoadVals(axes, L, q, load_type, a=a, b=b)
#

