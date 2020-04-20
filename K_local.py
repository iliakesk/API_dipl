#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 23:50:11 2020

@author: ilias
"""

#
#    def groupVals(self,group):
#        Group=self.Groups[group]
#        self.E=Group[0][0]
#        self.A=Group[0][1]
#        self.v=Group[0][2]
#        self.G=Group[0][3]
#        if self.G==0:
#            self.G = self.E/(2+2*self.v)
#        self.Ix=Group[0][4]
#        self.Iy=Group[0][5]
#        self.Iz=Group[0][6]
#        self.selfload=Group[0][7]
#        self.sectiontype = Group[1]
#        self.group=group
#
#
#30000000.0, 0.125, 0.0, 12500000.0, 0.00179, 0.0026, 0.000651, 0.0


model = {"specs":
            {"E":30000000, 
             "A":0.125, 
             "v":0, 
             "G":12500000, 
             "Ix":0.00179, 
             "Iy":0.0026, 
             "Iz":0.000651},        
        "nodes":
            {"1":{"x":0, "y":0, "z":0}, 
             "2":{"x":0, "y":0, "z":4}, 
             "3":{"x":0, "y":4, "z":4}, 
             "4":{"x":-4, "y":4, "z":8}},
        "members":
            {"1":{"snode":1, "enode":2},
             "2":{"snode":2, "enode":3},
             "3":{"snode":3, "enode":4}}
        }
print(model.nodes)        

#
#def Klocal():
#
#
#
#
#
#    
#    
#def beamLength(self, start, stop):
#    length = np.sqrt([(stop[0]-start[0])**2 + (stop[1]-start[1])**2 + (stop[2]-start[2])**2])
#    return length[0]
#
#
#K =np.zeros(shape=(12,12))
#        K[0][0]=K[6][6] = E*A/L
#        K[1][1]=K[7][7] = 12*E*Iz/(L**3)
#        K[2][2]=K[8][8] = 12*E*Iy/(L**3)
#        K[3][3]=K[9][9] = G*Ix/L
#        K[4][4]=K[10][10] = 4*E*Iy/L
#        K[5][5]=K[11][11] = 4*E*Iz/L
#        K[1][5]=K[5][1]=K[1][11]=K[11][1] = 6*E*Iz/(L**2)
#        K[7][11]=K[11][7]=K[5][7]=K[7][5] = -6*E*Iz/(L**2)
#        K[4][8]=K[8][4]=K[8][10]=K[10][8] = 6*E*Iy/(L**2)
#        K[2][4]=K[4][2]=K[2][10]=K[10][2] = -6*E*Iy/(L**2)
#        K[0][6]=K[6][0]=-K[0][0]
#        K[1][7]=K[7][1]=-K[1][1]
#        K[2][8]=K[8][2]=-K[2][2]
#        K[3][9]=K[9][3]=-K[3][3]
#        K[4][10]=K[10][4]=K[4][4]/2
#        K[5][11]=K[11][5]=K[5][5]/2
#        self.Klocal = K         