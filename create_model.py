#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 23:50:11 2020

@author: ilias
"""
import numpy as np
import element as el
import loads as lds
#import json_model
#
#model = json_model.model

#to object model tha pigainei sto class stucture kai apo ekei tha dinei mono oti
#xreiazetai sto class element
model = {
        "structure_type":"3dframe", #[2dframe, 3dframe, 2dtruss, 3dtruss, ........]
        "last_node_name":4,
        "last_element_name":3,
        "sections":
            {"1":{"name":"custom",
                  "specs":{ "E":30000000, 
                            "A":0.125, 
                            "v":0, 
                            "G":12500000, 
                            "Ix":0.00179, 
                            "Iy":0.0026, 
                            "Iz":0.000651}}},        
        "nodes":
            {"1":{"x":0, "y":0, "z":0, "load_ids":[]}, 
             "2":{"x":0, "y":0, "z":4, "load_ids":[1]}, 
             "3":{"x":0, "y":4, "z":4, "load_ids":[2]}, 
             "4":{"x":-4, "y":4, "z":8, "load_ids":[]}},
        "members":
            {"1":{"snode":1, "enode":2, "section":1, "theta":0, "edge_release":[0,0,0,0,0,0,0,0,0,0,0,0], "type":"generic"},
             "2":{"snode":2, "enode":3, "section":1, "theta":2.355, "edge_release":[0,0,0,0,0,0,0,0,0,0,0,0], "type":"generic"},
             "3":{"snode":3, "enode":4, "section":1, "theta":1.57, "edge_release":[0,0,0,0,0,0,0,0,0,0,0,0], "type":"generic"}},
        "supports":
            {"1":"111111"},                  
        "loads":
            {"1":{"load_type":"nodal",
                  "element":"2",
                  "value":[-20.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                  "group":"G"},
            "2":{"load_type":"nodal",
                  "element":"3",
                  "value":[0.0, 0.0, -20.0, 0.0, 0.0, 0.0],
                  "group":"G"},
#            "3":{"load_type":"E",#or T or R or Po etc
#                 "start":0,
#                 "stop":3,
#                 "element":"1",
#                 "value":[0.0, 0.0, -20.0, 0.0, 0.0, 0.0],
#                 "group":"G"}
            }
        }


members = {}
nodes = {}
supports = {}
nodal_load_sum = {}
member_loads = {}
#rotated_supports = {}
displacements = {}
#elastic_node = {}
#rigid_node = {}
#sections = {}



#populate the members dictionary
for key in model["members"].keys():
    name = key
    snode_name = model["members"][name]["snode"]
    enode_name = model["members"][name]["enode"]
    section_name = model["members"][name]["section"]
    
    snode = model["nodes"][str(snode_name)]
    enode = model["nodes"][str(enode_name)]
    section = model["sections"][str(section_name)]
    theta = model["members"][name]["theta"]
    releases = model["members"][name]["edge_release"]
    mtype = model["members"][name]["type"]
    
    members[name] = el.Element(name, snode_name, enode_name, snode, enode, section, theta, releases, mtype)

#populate the nodes and create the nodal_load_sum dictionary
for key in model["nodes"].keys():
    nodes[key] = model["nodes"][key]
    nodal_load_sum[key] = np.zeros(6)
    displacements[key] = np.zeros(6)



#populate the nodal_load_sum dictionary
for key in model["loads"].keys():
    if model["loads"][key]["load_type"] == "nodal":
        node = model["loads"][key]["element"]
        nodal_load_sum[node] += model["loads"][key]["value"]
    else:
        load_type = model["loads"][key]["load_type"]
        a = model["loads"][key]["start"]
        b = model["loads"][key]["stop"]
        q = model["loads"][key]["value"]
        member = model["loads"][key]["element"]
        L = members[member].L
        axes = members[member].axes
        snode_name = members[member].snode_name
        enode_name = members[member].enode_name
        
        lstart, lstop = lds.getEdgeValues(axes, L, q, load_type, a, b)
        #na dw ti xreiazetai  na ftiaxw gia na exw kai ena dictionary gia ta 
        #member loads gia na mporw na vrw meta ta internal forces
        nodal_load_sum[snode_name] -= lstart
        nodal_load_sum[enode_name] -= lstop
        
#populate the supports dictionary
for key in model["supports"].keys():
    supports[key] = model["supports"][key]





#calculate disps and reacts demo
def dispsReacts():
    nod_sum = len(nodes)
    s=6 # size 'supported' for the analysis of the matrices
    f = nod_sum*6-s #size 'free' for the analysis of the matrices
    K = K_global()
    tM = transformMatrix()
    Fall = forceMatrix()
    Dmatrix = displacementMatrix()
    Kt = np.dot(tM,np.dot(K,np.linalg.inv(tM)))#transformed K
    Fall_t = np.dot(tM,Fall)
    Dt = np.dot(tM,Dmatrix)
    Kff, Kfs, Ksf, Kss = Kt[:f,:f], Kt[:f,f:], Kt[f:,:f], Kt[f:,f:]#sliced Kff, Kfs etc
    Ff  = Fall_t[:f,]
    Fs1 = Fall_t[f:,]
    Ds = Dt[f:]
    try:
        Kffinv = np.linalg.inv(Kff)
    except:
        Kffinv = np.linalg.pinv(Kff)
        print('! ! ! ! ! In forcematrix of analysis3dframe 251 the pinv was used instead of inv')
    Df = np.dot(Kffinv,Ff-np.dot(Kfs,Ds))
    Fs = np.dot(Ksf,Df) + np.dot(Kss,Ds)

    Fs = Fs - Fs1
    print(Df)
    print(Fs)
#    DISPS = returnDisplacements(Df)
#    REACTS = returnReactions(Fs)
#    if thereΙsElasticNode():
#        DISPS,REACTS = calculateElasticNode(DISPS, REACTS)
#    return DISPS,REACTS


def K_global():
    '''Returns the stiffness matrix of the structure'''
#    rotatedsupports = st.Beam.RotatedSupports
#    elasticnode = st.Beam.ElasticNode
    nod_sum = len(nodes)*6
    Kall = np.zeros(shape=(nod_sum,nod_sum))
    
    for e in members.values():
        start = int(e.snode_name)-1
        stop = int(e.enode_name)-1
        Kel = e.K
        for i in range(6):
            for j in range(6):
                Kall[i+start*6][j+start*6] += Kel[i][j]
                Kall[i+start*6][j+stop*6] += Kel[i][j+6]
                Kall[i+stop*6][j+start*6] += Kel[i+6][j]
                Kall[i+stop*6][j+stop*6] += Kel[i+6][j+6]
#    if len(rotatedsupports)>0:
#        for node, angles in rotatedsupports.items():
#            rotv = supportRotationMatrix(node,angles)
#            Kall = np.dot(np.linalg.inv(rotv),np.dot(Kall,rotv))
#    if len(elasticnode)>0:
#        allnodes = sorted(nodes.keys(), key = lambda x: int(x[1:]))
#        for node in allnodes:
#            index = 6*(int(node[1:])-1)
#            if node in elasticnode.keys():
#                transfmatrix = st.elasticNode(node,elasticnode[node])
#                for l in range(6):
#                    Kall[index+l][index+l] += transfmatrix[l][l]
    return Kall




def displacementMatrix():
    n = len(nodes)*6
    Dmatrix = np.zeros(shape=(n,1))
    for node in sorted(nodes.keys(), key = lambda x: int(x)):
        if node in displacements.keys():
            index = 6*(int(node)-1)
            for i in range(6): 
                Dmatrix[index+i] = displacements[node][i]
    return Dmatrix




def transformMatrix():
    '''Returns the transformation matrix given the number of freedoms per node in general'''
    n = len(nodes)*6
    T = np.identity(n)
    transfMatrix = np.zeros(shape=(n,n))
    keynodes = sorted(supports.keys(), key = lambda x: int(x))
    numfreedom = []
    for i, node in enumerate(keynodes):
        for position, support in enumerate(supports[node]):
            if support=='1':
                numfreedom.append(position+6*i+1)
    count = 0            
    firstpos = n - len(numfreedom)
    for f in range(n):
        if f+1 not in numfreedom:
            transfMatrix[f-count]=T[f]
        else:
            transfMatrix[firstpos+count] = T[f]
            count+=1
    return transfMatrix



def forceMatrix():
    '''Return the force matrix for each group of loads. Need to be called for
    each group specifically (the load combinations will be formed after this)'''
 
    n = len(nodes)*6
    Pmatrix = np.zeros(shape=(n,1))
    P = nodal_load_sum

    i=0
    for node in sorted(P.keys(), key=lambda x: int(x)):
        Pmatrix[i] = P[node][0]
        Pmatrix[i+1] = P[node][1]
        Pmatrix[i+2] = P[node][2]
        Pmatrix[i+3] = P[node][3]
        Pmatrix[i+4] = P[node][4]
        Pmatrix[i+5] = P[node][5]
        i+=6
    return Pmatrix



dispsReacts()


