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
from enum import Enum
from scipy.optimize import linear_sum_assignment
from copy import  deepcopy

class ALGO_TYPE(Enum):
    NAIVE = 1
    HOP_DIVIDED = 2
    EPSILON_TRICK = 3
    BENCHMARK_1 = 4
    BENCHMARK_2= 5

class Multihop_Matching( Matching ):
    def __init__(self, W, delta, topo_file, traffic_file, routing_file = None ):
        super().__init__(topo_file, traffic_file, routing_file)

        self.W = W
        self.delta = delta
        self.edge_weights = None
        self.cur_matching = None
        self.cur_alpha = -1
        self.matching_history = []
        #self.packet_forwarding_history = []

        self.algo_type = ALGO_TYPE
        #self.cur_algo_type = ALGO_TYPE.NAIVE
        return

    def calculate_edge_weights(self):
        '''

        :return:
        '''
        self.edge_weights = np.zeros((self.n, self.n), dtype=np.float)
        for cur_node, flows in enumerate( self.current_next_hop_traffic):
            for cur_flow in flows:
                src, dest, next_hop, flow_val = cur_flow
                if self.topology[cur_node, next_hop]:
                    self.edge_weights[cur_node, next_hop] += flow_val #TODO: divide by path length
        return

    def find_M_alpha(self, max_duration = -1, hop_trick = False, epsilon_trick = False):
        '''

        :param max_duration:
        :return:
        '''
        best_score = -1
        best_matching = None
        best_alpha = -1
        alpha_set = list( set( self.edge_weights[ np.where( self.edge_weights > 0 )].flatten() ) )
        if len(alpha_set) <=0 : #no traffic left
            return best_matching, best_alpha

        if max_duration >= 0:
            alpha_set = [i  for i in alpha_set if i<=max_duration ]
            if len(alpha_set) == 0 and max_duration >0: #in case no alpha after above, just use the remaining time i,e max_duration
                alpha_set = [max_duration]

        for alpha in alpha_set:
            clipped_weights = np.clip( self.edge_weights , a_max=alpha, a_min=0 )
            row_indx, col_indx = linear_sum_assignment(- clipped_weights)
            matching_score = np.sum(clipped_weights[row_indx, col_indx]) / (1.0 * alpha + self.delta) #TODO: incorporate hop and e-trick
            if matching_score > best_score:
                best_score, best_alpha = matching_score, alpha
                best_matching = list( zip(row_indx, col_indx) )
        return best_matching, best_alpha


    def rank_flows(self, flow_list):
        '''

        :param flow_list:
        :return:
        '''
        return sorted( flow_list, key = lambda x: x[2], reverse=True )



    def route_flows(self):
        '''
        :return:
        '''
        snapshot_cur_matching = deepcopy(self.current_next_hop_traffic)
        for (n1, n2)  in self.cur_matching:
            if not self.topology[n1, n2]: #no connection, useless matching
                continue
            all_flows_between_n1_n2 = self.find_flows_between_nodes(n1, n2, nh_traffic= snapshot_cur_matching)
            ranked_flows_between_n1_n2 = self.rank_flows( all_flows_between_n1_n2 )
            remaining_time = self.cur_alpha
            for (src, dest, f) in ranked_flows_between_n1_n2:
                if remaining_time <= 0: break
                routable_flow =  min(f, remaining_time)
                remaining_time -= routable_flow
                self.forward_packets( cur_node=n1, src=src, dest=dest, in_flow_val= routable_flow )
        return

    def solve_multihop_routing_naive(self):
        '''

        :return:
        '''
        iteration_count = 0
        print("iteration: ", iteration_count, " a/t  ", self.cur_alpha, "/", self.cur_time)
        self.debug_pretty_print_init_traffic()
        while self.cur_time < self.W:
            iteration_count += 1
            remaining_time = self.W - self.cur_time
            self.calculate_edge_weights()
            m, alpha = self.find_M_alpha(max_duration=remaining_time)
            if m is None: break;
            self.cur_time += int( alpha )
            self.cur_matching, self.cur_alpha = m, alpha
            self.matching_history.append( (self.cur_matching, self.cur_alpha) )
            self.route_flows()

            print("iteration: ",iteration_count, "a/t ",self.cur_alpha,"/",self.cur_time)
            self.debug_pretty_print_current_traffic()
        return


    def solve_multihop_routing(self, algo_type):
        if algo_type == self.algo_type.NAIVE:
            self.solve_multihop_routing_naive()
        return

    def debug_pretty_print_current_traffic(self):
        '''
        :return:
        '''
        for val in self.cur_matching:
            n1, n2 = val
            print(n1, "-->",n2)
            for val in self.current_next_hop_traffic[n1]:
                src, dest, nh, f = val
                print("\t(", src,",", dest,")->",nh," f:",f)
        print("=======================================")
        return
    #-------end of class definition--------------#

if __name__ == '__main__':
    '''
    module test
    '''
