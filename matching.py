#!/usr/bin/env python
'''
Module: matching
List of entities:
Use:

'''
from __future__ import  print_function

__author__ = "Md Shaifur Rahman"
__email__ = "mdsrahman@cs.stonybrook.edu"
__date__ = "May 24, 2019"
__license__ = "GPL"
__copyright__ = "Copyright 2019, WINGS Lab, Stony Brook University"


#---start here--->
import numpy as np
from termcolor import colored
from scipy.optimize import linear_sum_assignment
from copy import  deepcopy

class Matching(object):
    def __init__(self, topo_file, traffic_file, routing_file = None):
        '''

        '''
        self.n = 0 #no of nodes/racks
        self.topology = None # connectivity matrix
        self.traffic = None #traffic demand matrix
        self.routing_table =  None #next hop matrix

        self.edge_weights = None #edge-weights calculated by the calculate_edge_weights(..) method
        self.current_traffic = None
        self.current_routing_table = None

        self.load_topology(topo_file)
        self.load_traffic(traffic_file)
        if routing_file is not None:
            self.load_routes(routing_file)

        return


    def load_topology(self, input_file):
        '''
        :param input_file: 
        :return: 
        '''
        is_bidirectional = True
        with open(input_file, 'r') as f:
            words = f.readline().split()
            self.n =  int( words[0] )
            if len(words) > 1:
                if 'c' == words[1].lower():
                    self.topology = np.ones((self.n, self.n), dtype=np.bool)
                    return
                elif 'u' == words[1].lower():
                    is_bidirectional = False
                elif 'b' == words[1].lower():
                    is_bidirectional =True
            #else
            self.topology = np.zeros((self.n, self.n), dtype=np.bool)
            for line in f:
                words = line.split()
                if len(words) < 2: continue #handle the blank space
                try:
                    i, j = int( words[0] )-1, int( words[1] )-1
                except:
                    print(colored( ("Error @load_topology(..): Invalid indices for i,j "), "red") )
                    exit(1)
                if i < self.n and i >= 0 and j < self.n and j >= 0:
                    self.topology[i, j] = True
                    if is_bidirectional:
                        self.topology[j, i] = True
        return

    def load_traffic(self, input_file, delimiter = ','):
        '''
        :param input_file:
        :return:
        '''
        try:
            self.traffic = np.loadtxt(input_file, delimiter=delimiter, dtype=np.longlong)
        except:
            print( colored( ("Error @load_traffic(..) Error in loading traffic matrix from file ", input_file), "red") )
            exit(1)
        if self.traffic.shape[0] > self.n or self.traffic.shape[1] > self.n:
            print( colored( ("Error @load_traffic(..) Traffic-matrix shape from file ", input_file,", does not match #-of nodes from topology file !!"), "red") )
            exit(1)
        self.current_traffic = deepcopy(self.traffic)
        return

    def load_routes(self, input_file, inter_node_delim = ',', src_dest_path_delim = ':'):
        '''

        :param input_file:
        :return:
        '''
        self.routing_table = [ [ [] for i in range(self.n) ] for j in range(self.n) ]

        with open(input_file, 'r') as f:
            for line in f:
                if line.isspace(): continue

                src_dest, intermediates = line.split(src_dest_path_delim)

                if src_dest.isspace() or  intermediates.isspace(): continue

                src_dest = src_dest.split( inter_node_delim )
                intermediates = intermediates.split( inter_node_delim )

                try:
                    i , j = int( src_dest[0] )-1, int( src_dest[1] )-1
                    if 0 <= i and i < self.n and 0 <=j and j < self.n:
                        next_hops = []
                        for k in intermediates:
                            k = int( k ) - 1
                            if 0 <= k and k < self.n  :
                                next_hops.append(k)
                            else:
                                raise Exception('')
                        self.routing_table[i][j] = next_hops
                except:
                    print( colored( ("Error @load_routes(..) Error in loading route matrix from file ", input_file), "red") )
        self.current_routing_table = deepcopy(self.routing_table)
        return

    def update_traffic(self, n1, n2,  ):

        return

    def calculate_edge_weights(self):
        '''
        abstract class
        :return:
        '''
        self.edge_weights = np.zeros( (self.n, self.n), dtype=np.float )

        return

    def get_bipartite_matching(self):
        '''
        :return:
        '''
        row_indx, col_indx = linear_sum_assignment(- self.edge_weights)

        return row_indx, col_indx

    #-------end of class definition--------------#

if __name__ == '__main__':
    '''
    module test
    '''
    topo_file = './data/synthetic/topology.txt'
    traffic_file = './data/synthetic/traffic.txt'
    routing_file = './data/synthetic/routing.txt'

    m = Matching(topo_file=topo_file, traffic_file=traffic_file, routing_file=routing_file)

    m.calculate_edge_weights()

    r, c = m.get_bipartite_matching()
    print(r, c)
