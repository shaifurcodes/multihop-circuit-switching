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


class Matching(object):
    def __init__(self, topo_file, traffic_file, routing_file, algo_type ):
        '''

        '''
        self.cur_algo_type = algo_type
        self.n = 0 #no of nodes/racks
        self.topology = None # connectivity matrix
        self.traffic = None #traffic demand matrix
        self.routing_table =  None #next hop matrix

        self.edge_weights = None #edge-weights calculated by the calculate_edge_weights(..) method
        self.current_next_hop_traffic = None #
            # 2D array each containing list of tuples (src, dest, next_node, flow_value, path_pointer, traversed_hops, remaining_hops, chosen_path_index)
        self.stat_routed_flows = None

        self.cur_time = 0
        self.max_hop = 0
        self.multipath_factor = 1
        self.debug_overflow_amount = 0

        self.load_topology(topo_file)
        self.load_traffic(traffic_file)

        if routing_file is not None:
            self.load_routes(routing_file)
            self.init_current_next_hop_traffic()
        else:
            exit(1)


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
                    np.fill_diagonal(self.topology, False) #avoid loops
                    return
                elif 'u' == words[1].lower():
                    is_bidirectional = False
                elif 'b' == words[1].lower():
                    is_bidirectional =True

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

        np.fill_diagonal(self.topology, False) #avoid loops
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
        self.max_hop = 1
        path_count = np.zeros((self.n, self.n), dtype=np.int)
        with open(input_file, 'r') as f:
            for line in f:
                if line.isspace(): continue
                src_dest, intermediates = line.split(src_dest_path_delim)
                if src_dest.isspace(): continue
                try:
                    src_dest = src_dest.split(inter_node_delim)
                    i , j = int( src_dest[0] ), int( src_dest[1] )
                    #----count multipath----#
                    path_count[i, j] +=  1
                    #----------------------#
                    if intermediates.isspace():
                        self.routing_table[i][j].append([])  # add empty path
                        continue
                    intermediates = intermediates.split(inter_node_delim)
                    if 0 <= i and i < self.n and 0 <=j and j < self.n:
                        next_hops = []
                        for k in intermediates:
                            k = int( k )
                            if 0 <= k and k < self.n  :
                                next_hops.append(k)
                            else:
                                raise Exception('')
                        self.routing_table[i][j].append(next_hops)
                        self.max_hop = max(self.max_hop, 1+len(next_hops))
                except:
                    print( colored( ("Error @load_routes(..) Error in loading route matrix from file ", input_file), "red") )
                    exit(1)
        self.multipath_factor = np.max(path_count)

        if self.cur_algo_type == 'upperbound':
            #----------modify routing table if upperbound, convert all routes to single hop------#
            #------first add the relevant flows:
            upperbound_traffic = np.zeros_like(self.traffic, dtype=np.longlong)
            for i in range(self.n):
                for j in range(self.n):
                    if self.traffic[i, j] >0:
                        cur_path  = [i]+ self.routing_table[i][j][0] + [j] #add src and dest
                        for k in range( len(cur_path)-1):
                            x, y = cur_path[k], cur_path[k+1]
                            upperbound_traffic[x, y] = self.traffic[i, j]

            #-----now fix the upperbound routes all should be empty []:---
            upperbound_routing_table = [[[] for i in range(self.n)] for j in range(self.n)]

            #----now replace the traffic and routes for upperbound---#
            self.traffic = upperbound_traffic
            self.routing_table = upperbound_routing_table
        # if self.cur_algo_type == 'benchmark_1':
        #     #split the traffic per hop
        #     for i in
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
                    for cur_path_indx, cur_path in enumerate(self.routing_table[i][j]):
                        cur_path_length = len(cur_path)
                        if len(cur_path) > 0:
                            next_hop = cur_path[0]
                        else:
                            next_hop = j
                        traversed_hops = 0
                        remaining_hops = 1+cur_path_length
                        self.current_next_hop_traffic[i].append( (i, j, next_hop, flow_val, traversed_hops, remaining_hops,cur_path_indx ) )
        return

    def debug_pretty_print_init_traffic(self):
        '''

        :return:
        '''
        for n1 in range(self.n):
            print("node: ",n1)
            for val in self.current_next_hop_traffic[n1]:
                src, dest, nh, f, _, _ , p_indx= val
                print("\t(", src,",", dest,")->",nh," f:",f, " pid: ",p_indx)
        print("=======================================")
        return

    def find_next_hop(self, src, dest, cur_hop, path_indx):
        '''

        :param src:
        :param dest:
        :param cur_hop:
        :return:
        '''
        if self.routing_table is None:
            print("routing table still empty!!")
            exit(1)
        else:
            routing_path = self.routing_table[src][dest][path_indx]
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

    def find_flows_between_nodes(self, cur_node, next_node, nh_traffic = None):
        '''

        :param cur_node:
        :param next_node:
        :param nh_traffic:
        :return:
        '''
        cur_nh_traffic = nh_traffic
        if cur_nh_traffic is None:
            cur_nh_traffic = self.current_next_hop_traffic
        flow_list = []
        for val in  cur_nh_traffic[cur_node] :
            i , j , nh, f, t_thop, r_hop, path_indx = val
            if nh == next_node:
                flow_list.append((i, j, nh, f, t_thop, r_hop, path_indx))
        return flow_list

    def find_flow_indx_at_node(self, cur_node, src, dest, path_indx, nh_traffic = None):
        '''

        :param cur_node:
        :param src:
        :param dest:
        :param nh_traffic:
        :return:
        '''
        src_dest_indx  = -1
        cur_nh_traffic = nh_traffic
        if cur_nh_traffic is None:
            cur_nh_traffic = self.current_next_hop_traffic
        for indx, val in  enumerate( cur_nh_traffic[cur_node] ) :
            i , j , nh, f, _, _, p_indx = val
            if i==src and j==dest and path_indx == p_indx:
                src_dest_indx = indx
                break
        return src_dest_indx

    def remove_flows(self, cur_node, src, dest, f_reduce, path_indx):
        '''

        :param cur_node:
        :param src:
        :param dest:
        :param f_reduce:
        :return:
        '''
        new_flow_list = []
        for cur_indx, val in enumerate( self.current_next_hop_traffic[cur_node] ):
            i , j , nh, f, t_hop, r_hop, p_indx = val
            if i== src and j==dest and path_indx == p_indx:
                if f > f_reduce:
                    new_flow_list.append( (i , j , nh, f - f_reduce, t_hop, r_hop, p_indx) )
            else:
                new_flow_list.append(val)
        self.current_next_hop_traffic[cur_node]  = new_flow_list
        return

    def remove_flows_at_src(self, src, dest, except_p_indx, save_one_hop):
        '''

        :param src:
        :param dest:
        :param except_p_indx:
        :param save_one_hop:
        :return:
        '''
        new_flow_list = []
        for cur_indx, val in enumerate( self.current_next_hop_traffic[src] ):
            i , j , nh, f, t_hop, r_hop, p_indx = val
            if i== src and j==dest and p_indx != except_p_indx:
                if save_one_hop:
                    if len(self.routing_table[src][dest][p_indx])<=0: #1-hop flow
                        new_flow_list.append( val )
            else:
                new_flow_list.append(val)
        self.current_next_hop_traffic[src]  = new_flow_list
        return

    def find_all_intermediate_nodes(self, src, dest, path_indx = -1):
        '''

        :param src:
        :param dest:
        :param path_indx:
        :return:
        '''
        node_list = []
        if path_indx >=0:
            node_list = self.routing_table[src][dest][path_indx]
        else:
            for path in self.routing_table[src][dest]:
                node_list.extend(path)
            node_list = list(set(node_list))
        return  node_list

    def forward_packets(self, cur_node, src, dest, path_indx, in_flow_val=-1, delete_duplicate_flows = 0):
        '''
        route routed_flow_val amount of packets from n to the next hop for (src, dest) flow
        if routed_flow_val is -1, then flow whatever flow left at that node
        :param cur_node: current node
        :param src:
        :param dest:
        :param in_flow_val:
        :return:
        '''

        cur_flow_val, cur_next_hop, cur_indx, cur_traversed_hop_count, cur_remaining_hop_count, cur_p_indx =\
                                -1, -1, -1, -1, -1, -1
        for indx, val in enumerate( self.current_next_hop_traffic[cur_node] ):
            i , j , nh, f, traversed_hop_count, remaining_hop_count, p_indx = val
            if i== src and j==dest and path_indx == p_indx:
                cur_indx = indx
                cur_flow_val, cur_next_hop,  cur_traversed_hop_count, cur_remaining_hop_count, cur_p_indx  = \
                           f,           nh,      traversed_hop_count,     remaining_hop_count, p_indx
                break

        if cur_flow_val == -1:
            print(colored(("Warning @route_traffic(..) current (src, dest):", src, dest," not found in cur_node: ", cur_node), "red"))
            return
        #if no link, ignore
        if not self.topology[cur_node, cur_next_hop]:
            return
        #else---
        routed_flow_val = in_flow_val

        if routed_flow_val ==-1:
            routed_flow_val = cur_flow_val


        if dest == cur_next_hop: #destination reached for this flow
            self.collect_stat(src, dest, routed_flow_val)


        else: #--update entry for the next hop_node
            src_dest_indx = self.find_flow_indx_at_node(cur_node=cur_next_hop,
                                                        src=src,
                                                        dest=dest,
                                                        path_indx=cur_p_indx)
            if src_dest_indx >=0:
                src, dest, nh, f,  traversed_hop_count, remaining_hop_count, p_indx = \
                    self.current_next_hop_traffic[cur_next_hop][src_dest_indx]
                self.current_next_hop_traffic[cur_next_hop][src_dest_indx] = \
                    (src, dest, nh, f+routed_flow_val, traversed_hop_count, remaining_hop_count, p_indx)
            else:
                new_next_hop = self.find_next_hop(src, dest, cur_next_hop, cur_p_indx )
                self.current_next_hop_traffic[cur_next_hop].append( (src, dest, new_next_hop,
                                                                     routed_flow_val,
                                                                     cur_traversed_hop_count+1,
                                                                     cur_remaining_hop_count-1,
                                                                     cur_p_indx) )
        #update entry for the current hop
        if cur_flow_val - routed_flow_val >0 : #flow remaining
            self.current_next_hop_traffic[cur_node][cur_indx] =\
                src, dest, cur_next_hop, (cur_flow_val - routed_flow_val), cur_traversed_hop_count, \
                cur_remaining_hop_count, cur_p_indx
        else:
            del self.current_next_hop_traffic[cur_node][cur_indx]

        if cur_node==src:
            if delete_duplicate_flows == 1:
                self.remove_flows_at_src(src, dest, cur_p_indx, False)
            elif delete_duplicate_flows == 2 :
                self.remove_flows_at_src(src, dest, cur_p_indx, True)

        if cur_next_hop==dest and delete_duplicate_flows == 2:
                self.remove_flows_at_src(src, dest, cur_p_indx, True)
        return

    def collect_stat(self, src, dest, flow_val):
        '''

        :param src:
        :param dest:
        :param flow:
        :return:
        '''

        self.stat_routed_flows[src, dest] += flow_val
        print("DEBUG: flow completed (src, dest, flow): ",src,dest, flow_val)
        if self.stat_routed_flows[src, dest] > self.traffic[src, dest]:
            #---there are cases with backtracking where duplicates flow packets exist
            #----i.e duplicate flows reach destination in the same iteration, so ignore  them
            self.debug_overflow_amount += self.stat_routed_flows[src, dest]- self.traffic[src, dest]
            self.stat_routed_flows[src, dest] =  self.traffic[src, dest]
        return

    def find_demand_met(self):
        '''

        :return:
        '''
        total_demand = np.sum(self.traffic)
        demand_met = np.sum(self.stat_routed_flows)
        return 1.*demand_met/total_demand
    #-------end of class definition--------------#

if __name__ == '__main__':
    '''
    module test
    '''
    base_file_name = './data/synthetic/multipath_1'
    topo_file = base_file_name+'.topology.txt'
    traffic_file = base_file_name+'.traffic.txt'
    routing_file = base_file_name+'.routing.txt'
    m = Matching(topo_file=topo_file, traffic_file=traffic_file, routing_file=routing_file)

