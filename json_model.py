#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 22:35:08 2020

@author: ilias
"""

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
