#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 18:38:48 2020

@author: ilias
"""

import numpy

#to class Element tha pairnei mono oti xreiazetai diladi ta tou kathe stoixeioy
class Element():
    def __init__(self, name, snode_name, enode_name, snode, enode, section, theta, releases,mtype="generic"):
#        print(snode)
        self.name = name
        self.snode_name = int(snode_name)
        self.enode_name = int(enode_name)
        self.snode = numpy.array([snode["x"], snode["y"], snode["z"]])
        self.enode = numpy.array([enode["x"], enode["y"], enode["z"]])
        self.E = section["specs"]["E"]
        self.A = section["specs"]["A"]
        self.v = section["specs"]["v"]
        self.G = section["specs"]["G"]
        self.Ix = section["specs"]["Ix"]
        self.Iy = section["specs"]["Iy"]
        self.Iz = section["specs"]["Iz"]
        self.theta = theta
        self.releases = releases
        self.mtype = mtype #member type for the design (beam, column etc) -to review
        
        self.L = self.length(self.snode, self.enode)
        self.axes = self.localAxes(self.snode, self.enode, self.theta)
        self.rotv, self.rotMatrix = self.rotationMatrix(self.axes)
        
        self.getKEl(self.L, self.A, self.E, self.G, self.Ix, self.Iy, self.Iz,self.v)
        


    def length(self, snode, enode):
        l = numpy.sqrt([(enode[0]-snode[0])**2 + (enode[1]-snode[1])**2 + (enode[2]-snode[2])**2])
        return l[0]
    
    
    
    def rotationMatrix(self, axes):
        '''Returns the rotation matrix of the element with axes being a tuple
        axes = (x_local, y_local, z_local)'''
        rotv = numpy.zeros(shape=(3,3))
        rotv[0],rotv[1],rotv[2]=axes[0],axes[1],axes[2]
        for l,line in enumerate(rotv):
            for v,value in enumerate(line):
                if abs(value)<1e-8:
                    rotv[l][v]=0
        rotationMatrix = numpy.zeros(shape=(12,12))
        for count in [0,3,6,9]:
            for i in range(3):
                for j in range(3):
                    rotationMatrix[i+count][j+count] = rotv[i][j]
        return rotv, rotationMatrix
    
    
    def rot_theta(self, axis, theta):
        axis = axis/numpy.sqrt(numpy.dot(axis, axis))
        a = numpy.cos(theta/2)
        b, c, d = -axis*numpy.sin(theta/2)
        return numpy.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                      [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                      [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])
    
    
    def localAxes(self, snode, enode, theta):
        '''Returns the unit vectors of the local axes of the member'''
        u = numpy.array([enode[i]-snode[i] for i in range(3)])
        x_local = u/numpy.linalg.norm(u)
        reflect = numpy.array([[-1,0,0],[0,-1,0],[0,0,1]])#<---
        uz = numpy.dot(u,reflect)
        if u[2]>0:
            if u[0]==0 and u[1]==0:
                z_local = numpy.array([-1,0,0])
                y_local = numpy.cross(z_local,u)/numpy.linalg.norm(numpy.cross(u,z_local))
            else:
                y_local = numpy.cross(uz,u)/numpy.linalg.norm(numpy.cross(u,uz))
                z_local = numpy.cross(u,y_local)/numpy.linalg.norm(numpy.cross(y_local,u))
        elif u[2]<0:
            if u[0]==0 and u[1]==0:
                z_local = numpy.array([1,0,0])
                y_local = numpy.cross(z_local,u)/numpy.linalg.norm(numpy.cross(u,z_local))
            else:
                y_local = numpy.cross(u,uz)/numpy.linalg.norm(numpy.cross(uz,u))
                z_local = numpy.cross(u,y_local)/numpy.linalg.norm(numpy.cross(y_local,u))
        elif u[2]==0:
            z_local = numpy.array([0,0,1])
            y_local = numpy.cross(z_local,u)/numpy.linalg.norm(numpy.cross(u,z_local))
        if theta!=0:
            rot_x = self.rot_theta(x_local, theta)
            y_local = numpy.dot(y_local, rot_x)
            z_local = numpy.dot(z_local, rot_x)
        return x_local, y_local, z_local
    
    
    
    def getKEl(self, L, A, E, G, Ix, Iy, Iz,v):        
        '''Returns the stiffness matrix of the member'''
        if G==None or G==0:
            G=E/(2+2*v)
        R = self.rotMatrix
        RT = numpy.linalg.inv(R)
        K = numpy.zeros(shape=(12,12))
        K[0][0]=K[6][6] = E*A/L
        K[1][1]=K[7][7] = 12*E*Iz/(L**3)
        K[2][2]=K[8][8] = 12*E*Iy/(L**3)
        K[3][3]=K[9][9] = G*Ix/L
        K[4][4]=K[10][10] = 4*E*Iy/L
        K[5][5]=K[11][11] = 4*E*Iz/L
        K[1][5]=K[5][1]=K[1][11]=K[11][1] = 6*E*Iz/(L**2)
        K[7][11]=K[11][7]=K[5][7]=K[7][5] = -6*E*Iz/(L**2)
        K[4][8]=K[8][4]=K[8][10]=K[10][8] = 6*E*Iy/(L**2)
        K[2][4]=K[4][2]=K[2][10]=K[10][2] = -6*E*Iy/(L**2)
        K[0][6]=K[6][0]=-K[0][0]
        K[1][7]=K[7][1]=-K[1][1]
        K[2][8]=K[8][2]=-K[2][2]
        K[3][9]=K[9][3]=-K[3][3]
        K[4][10]=K[10][4]=K[4][4]/2
        K[5][11]=K[11][5]=K[5][5]/2
        self.Klocal = K
        self.K = numpy.dot(RT,numpy.dot(K,R))    
    
    def getKElReleased(self, Rr, releases):
        '''Returns the stiffness matrix of the member if the member has any kind of edge release'''
        R = self.rotMatrix
        RT = numpy.linalg.inv(R)
        K = self.Klocal
        RrT = numpy.linalg.inv(Rr)
        Kr = numpy.dot(RrT,numpy.dot(K,Rr))
        Krglob = numpy.dot(RT,numpy.dot(Kr,R))
        self.rels, self.tM = transformReleases(releases)
        tMT = numpy.linalg.inv(self.tM)
        Kt = numpy.dot(self.tM,numpy.dot(Krglob,tMT))

        self.Kt = Kt
        tKt11 = numpy.linalg.pinv(Kt[:self.rels,:self.rels])
        Kt12 = Kt[:self.rels,self.rels:]
        Kt21 = Kt[self.rels:,:self.rels]
        Kt22 = Kt[self.rels:,self.rels:]
        Kr = Kt22 - numpy.dot(Kt21,numpy.dot(tKt11,Kt12))

        indexes = []
        for i, r in enumerate(releases):
            if r==1:
                indexes.append(i)
        for i in indexes:
            Kr = numpy.insert(Kr, i, 0, axis=1)
        for i in indexes:
            Kr = numpy.insert(Kr, i, numpy.zeros(12), axis=0)
#        self.Klocal = numpy.dot(Rr,numpy.dot(K,RrT))
        self.K = Kr


    def releaseEdge(self):
       #na elegjw an xreiazetai to angles. an xreiazetai tote prepei 
       #na diamorfwsw analoga to self.releases
        release = self.release
        release_angles = release[2]#[0,pi/3,0,0,0,0]
        releases = release[1]#[1,0,0,0,0,0]
        self.Rr = edgeRotationMatrix(release_angles)
        self.getKElReleased(self.Rr, releases)


    
    
#thelei diamorfwsi    
def transformReleases(releases):#releases =[0,0,0,0,0,0  ,  0,0,0,0,1,0]
    T = numpy.zeros(shape=(12,12))
    count = 0
    num_of_rel = 0
    for i, r in enumerate(releases):
        if r==1:
            T[count][i]=1
            count+=1
            num_of_rel+=1 #mipws thelei kapoio midenismo tou count ????????
    for i, r in enumerate(releases):
        if r==0:
            T[count][i]=1  #kati den paei kala
            count+=1
    return num_of_rel, T    
    
    
    
def edgeRotationMatrix(angles = (0,0,0)):
    phix = angles[0]                  
    phiy = angles[1]
    phiz = angles[2]
    rx = numpy.array([[1,0,0], [0, numpy.cos(phix), numpy.sin(phix)], [0, -numpy.sin(phix), numpy.cos(phix)]])
    ry = numpy.array([[numpy.cos(phiy), 0, numpy.sin(phiy)], [0,1,0], [-numpy.sin(phiy), 0, numpy.cos(phiy)]])
    rz = numpy.array([[numpy.cos(phiz), numpy.sin(phiz),0], [-numpy.sin(phiz), numpy.cos(phiz),0], [0,0,1]])
    rotv = numpy.dot(rx,numpy.dot(ry,rz))
    for l,line in enumerate(rotv):
        for v,value in enumerate(line):
            if abs(value)<1e-8:
                rotv[l][v]=0
    rotationMatrix = numpy.zeros(shape=(12,12))
    for count in [0,3,6,9]:
        for i in range(3):
            for j in range(3):
                rotationMatrix[i+count][j+count] = rotv[i][j]
    return rotationMatrix
