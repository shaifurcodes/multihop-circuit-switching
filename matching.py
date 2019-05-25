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
        self.current_next_hop_traffic = None # 2D array each containing list of tuples (src, dest, next_node, flow_value)
        self.stat_routed_flows = None

        self.alpha = 0

        self.load_topology(topo_file)
        self.load_traffic(traffic_file)
        if routing_file is not None:
            self.load_routes(routing_file)
            self.init_current_next_hop_traffic()
        else:
            self.init_c_n_h_traffic_without_input()

        self.stat_routed_flows = np.zeros((self.n, self.n), dtype=np.longlong)
        self.stat_fct = np.zeros((self.n, self.n), dtype=np.longlong)
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
                    exit(1)
        return

    def update_traffic(self, n1, n2,  ):

        return

    def init_current_next_hop_traffic(self):
        '''

        :return:
        '''
        self.current_next_hop_traffic = [ [] for i in range(self.n) ]
        for i in range(self.n):
            for j in range(self.n):
                flow_val = self.traffic[i, j]
                next_hop = -1
                if flow_val > 0:
                    if len(self.routing_table[i][j]) > 0:
                        next_hop = self.routing_table[i][j][0]
                    else:
                        next_hop = j
                    self.current_next_hop_traffic[i].append( (i, j, next_hop, flow_val) )

        return

    def init_c_n_h_traffic_without_input(self):
        '''

        :return:
        '''
        #TODO: without input file, generate intial routes/paths for current_next_hop_traffic matrix
        return

    def find_next_hop(self, src, dest, cur_hop):
        '''

        :param src:
        :param dest:
        :param cur_hop:
        :return:
        '''
        if self.routing_table is None:
            #TODO: handle missing routing table (generate route as well)
            print("routing table still empty!!")
            exit(1)
        else:
            routing_path = self.routing_table[src][dest]
            if len(routing_path) > 0: #indirect routing
                cur_hop_indx = -1
                try:
                    cur_hop_indx = routing_path.index(cur_hop)
                except  ValueError:
                    print (colored("Error @find_next_hop(..) next hop not found for (src,dest,curhop):",src,dest,cur_hop), "red")
                    exit(1)
                if cur_hop_indx+1 >= len(routing_path):
                    return dest
                else:
                    return routing_path[cur_hop_indx+1]
            else: #direct routing
                return dest
        return None # fix it for non-existent routing table

    def route_traffic(self, cur_node, src, dest, in_flow_val=-1):
        '''
        route routed_flow_val amount of packets from n to the next hop for (src, dest) flow
        if routed_flow_val is -1, then flow whatever flow left at that node
        :param cur_node: current node
        :param src:
        :param dest:
        :param in_flow_val:
        :return:
        '''
        cur_flow_val, cur_next_hop, cur_indx = -1, -1, -1
        for indx, val in enumerate( self.current_next_hop_traffic[cur_node] ):
            i , j , nh, f = val
            if i== src and j==dest:
                cur_flow_val, cur_next_hop, cur_indx = f, nh, indx
                break
        if cur_flow_val == -1:
            print(colored(("Error @route_traffic(..) current (src, dest):", src, dest," not found in cur_node: ", cur_node), "red"))
            return
        #else---
        routed_flow_val = in_flow_val
        if routed_flow_val ==-1:
            routed_flow_val = cur_flow_val


        if dest == cur_next_hop: #destination reached for this flow
            self.collect_stat(src, dest, routed_flow_val)
        else: #--update entry for the next hop
            new_next_hop = self.find_next_hop(src, dest, cur_next_hop)
            self.current_next_hop_traffic[cur_next_hop].append( (src, dest, new_next_hop, routed_flow_val) )

        remaining_flow = cur_flow_val - routed_flow_val

        #---updat entry for cur hop
        if remaining_flow>0: #update for remaining flows
            self.current_next_hop_traffic[cur_node][cur_indx] =  (src, dest, cur_next_hop, remaining_flow)
        else: #remove entry for current hop
            del self.current_next_hop_traffic[cur_node][cur_indx]
        return

    def collect_stat(self, src, dest, flow_val):
        '''

        :param src:
        :param dest:
        :param flow:
        :return:
        '''
        self.stat_routed_flows[src, dest] += flow_val
        if self.stat_routed_flows[src, dest] == self.traffic[src, dest]:
            self.stat_fct[src, dest] = self.alpha
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

    print(m.current_next_hop_traffic)

    m.alpha = 2
    m.route_traffic(0, 0, 1)
    m.route_traffic(0, 0, 3)
    m.route_traffic(1, 1, 2)
    m.route_traffic(1, 1, 3)
    m.alpha = 5
    print(m.current_next_hop_traffic)

    m.route_traffic(3, 0, 1)
    m.route_traffic(1, 0, 3)
    m.alpha = 10
    print(m.current_next_hop_traffic)

    m.route_traffic(2, 0, 1)
    m.alpha = 20
    print(m.current_next_hop_traffic)

    print(m.stat_routed_flows)
    print(m.stat_fct)
    #m.calculate_edge_weights()

    #r, c = m.get_bipartite_matching()
    #print(r, c)
