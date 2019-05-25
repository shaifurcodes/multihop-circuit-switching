#!/usr/bin/env python
'''
Module: multihop_matching
List of entities:
Use:

'''
from __future__ import  print_function
from builtins import  super

__author__ = "Md Shaifur Rahman"
__email__ = "mdsrahman@cs.stonybrook.edu"
__date__ = "May 24, 2019"
__license__ = "GPL"
__copyright__ = "Copyright 2019, WINGS Lab, Stony Brook University"


#---start here--->
from matching import  Matching
import  numpy as np

class Multihop_Matching( Matching ):
    def __init__(self, topo_file, traffic_file, routing_file = None ):
        super().__init__(topo_file, traffic_file, routing_file)
        return

    def calculate_edge_weights(self):
        super().calculate_edge_weights()
        #TODO: complete weight calculations
        return 
    #-------end of class definition--------------#
    
if __name__ == '__main__':
    '''
    module test
    '''
    topo_file = './data/synthetic/topology.txt'
    traffic_file = './data/synthetic/traffic.txt'
    routing_file = './data/synthetic/routing.txt'

    mmatching = Multihop_Matching(topo_file=topo_file, traffic_file=traffic_file, routing_file=routing_file)
    mmatching.calculate_edge_weights()
    r, c = mmatching.get_bipartite_matching()
    print (r, c)