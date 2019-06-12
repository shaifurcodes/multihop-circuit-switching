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
from timeit import default_timer as mytimer
from copy import deepcopy


class Multihop_Matching( Matching ):
    def __init__(self, W, delta, topo_file, traffic_file, algo_type, routing_file = None ):
        super().__init__(topo_file, traffic_file, routing_file)

        self.W = W
        self.delta = delta
        self.edge_weights = None
        self.cur_matching = None
        self.cur_alpha = -1
        self.matching_history = []
        #self.packet_forwarding_history = []

        self.cur_algo_type = algo_type
        self.forward_pkt_flag = 0
        if self.cur_algo_type == 'multipath':
            self.forward_pkt_flag = 1
        elif self.cur_algo_type == 'multipath_wbt':
            self.forward_pkt_flag = 2

        self.stat_matching_count = 0
        self.stat_matching_time = 0
        return

    def find_M_alpha(self, max_duration = -1 ):
        '''

        :param max_duration:
        :return:
        '''
        start_t2 = mytimer()
        best_score = -1
        best_matching = None
        best_alpha = -1
        alpha_set = self.generate_list_of_alphas()
        if len(alpha_set) <=0 : #no traffic left
            return best_matching, best_alpha

        if max_duration >= 0:
            alpha_set = [i  for i in alpha_set if i<=max_duration ]
            if len(alpha_set) == 0 and max_duration >0: #in case no alpha after above, just use the remaining time i,e max_duration
                alpha_set = [max_duration]

        for alpha in alpha_set:
            if alpha<=0:
                continue
            self.stat_matching_count+=1
            start_t = mytimer()

            start_t_w = mytimer()
            self.edge_weights = self.calculate_edge_weights(alpha)
            elapsed_t_w = mytimer() - start_t_w


            clipped_weights = np.clip( self.edge_weights , a_max=alpha, a_min=0 )

            start_t_m = mytimer()
            row_indx, col_indx = linear_sum_assignment(- clipped_weights)
            elapsed_t_m = mytimer() - start_t_m

            matching_score = np.sum(clipped_weights[row_indx, col_indx]) / (1.0 * alpha + self.delta)


            if matching_score > best_score:
                best_score, best_alpha = matching_score, alpha
                best_matching = list( zip(row_indx, col_indx) )
            self.stat_matching_time +=(mytimer() - start_t)
            #print("debug: w_time, m_time:", 1000*elapsed_t_w, 1000*elapsed_t_m)
        elapsed_t = mytimer() - start_t2
        print("debug: elapsed time:", elapsed_t, " seconds ", " #-alphas: ", len(alpha_set) )
        print("debug: list of alphas: ", alpha_set)
        return best_matching, best_alpha

    def virtual_route_flows(self):
        '''

        :return:
        '''
        # snapshot_cur_matching =  [ [] for i in range(self.n) ]
        # for i in range(self.n):
        #     for val in self.current_next_hop_traffic[i]:
        #         i, j, nh, f, path_ptr, t_hop, r_hop  =  val
        #         snapshot_cur_matching[i].append( (i, j, nh, f, path_ptr, t_hop, r_hop) )

        snapshot_cur_matching = deepcopy(self.current_next_hop_traffic)
        for (n1, n2)  in self.cur_matching:
            if not self.topology[n1, n2]: #no connection, useless matching
                continue
            all_flows_between_n1_n2 = self.find_flows_between_nodes(n1, n2, nh_traffic= snapshot_cur_matching)
            grouped_flow_list = self.group_flows_by_hops(all_flows_between_n1_n2)

            remaining_time = self.cur_alpha
            is_time_remaining = True
            for cur_grouped_flows in grouped_flow_list[1:]:
                for cur_flow in cur_grouped_flows:
                    if remaining_time <= 0:
                        is_time_remaining = False
                        break
                    src, dest, _, f, _, _  = cur_flow
                    routable_flow =  min(f, remaining_time)
                    remaining_time -= routable_flow

                    self.forward_packets( cur_node=n1, src=src, dest=dest, in_flow_val= routable_flow,
                                          delete_duplicate_flows=self.forward_pkt_flag)
                if not is_time_remaining:
                    break
        return

    def group_flows_by_hops(self, flow_list):
        '''

        :return:
        '''
        grouped_flow_list = [[] for i in range(self.max_hop+1)] #0-th index useless, accounting for one-to-one index mapping
        for val in flow_list:
            i, j, nh, f, t_hop, r_hop =  val
            grouped_flow_list[r_hop].append(val)
        return grouped_flow_list

    def find_edge_weight(self, n1, n2, alpha):
        '''

        :return:
        '''
        if not self.topology[n1, n2]:
            return 0
        flow_list = self.find_flows_between_nodes(n1, n2)
        if (len(flow_list)) == 0:
            return 0

        #else
        edge_benefit = 0
        sum_flows = 0
        grouped_flow_list = self.group_flows_by_hops(flow_list)
        alpha_exceeded =  False
        for hop_count, cur_group in enumerate( grouped_flow_list ):
            if hop_count <=0:
                continue
            cur_flows = sorted(cur_group, key = lambda x: x[3])
            for flow in cur_flows:
                _, _,_, flow_val,_,_ = flow
                routable_flows = 0
                if flow_val > alpha-sum_flows:
                    alpha_exceeded = True
                    routable_flows  = alpha - sum_flows
                else:
                    routable_flows = flow_val
                edge_benefit += 1.*routable_flows/hop_count
                if alpha_exceeded:
                    break #break for the innner for loop
            if alpha_exceeded:
                break #break for the outer for loop
        return edge_benefit

    def find_list_of_next_hops(self, n_indx):
        '''

        :param n_indx:
        :return:
        '''
        next_hops = []
        for val in self.current_next_hop_traffic[n_indx]:
            src, dest, nh, f, t_hop, r_hop, p_indx =  val
            next_hops.append(nh)
        return list(set(next_hops))

    def find_alpha_values_for_node(self, n_indx):
        '''

        :param n_indx:
        :return:
        '''
        list_of_next_hops = sorted( self.find_list_of_next_hops(n_indx) )

        sum_alpha_per_hop = [ [] for i in range(len(list_of_next_hops)) ]

        for val in sum_alpha_per_hop:
            for h in range(self.max_hop+1):
                val.append(0)


        for val in self.current_next_hop_traffic[n_indx]:
            i, j, nh, flow, t_hop, r_hop, p_indx = val
            nh_indx = list_of_next_hops.index(nh)
            if r_hop > 0:
                sum_alpha_per_hop[nh_indx][r_hop] += flow

        list_of_alphas = []
        for cur_nh in list_of_next_hops:
            cur_next_hop_alpha = [0 for i in range(self.max_hop+1)]

            for i in range(self.max_hop+1):
                for j in range(i, self.max_hop+1):
                    cur_nh_indx = list_of_next_hops.index(cur_nh)
                    cur_next_hop_alpha[j] += sum_alpha_per_hop[cur_nh_indx][i]
            list_of_alphas.extend(cur_next_hop_alpha)

        return  list_of_alphas

    def generate_list_of_alphas(self):
        '''

        :return:
        '''
        alpha_list = []
        for n_indx in range(self.n):
            n_alphas = self.find_alpha_values_for_node(n_indx)
            n_alphas = set(n_alphas)
            if 0 in n_alphas:
                n_alphas.remove(0)
            alpha_list.extend( list(n_alphas) )
        alpha_list = list( set(alpha_list ))
        return alpha_list

    def calculate_edge_weights(self, alpha):
        '''

        :return:
        '''
        edge_weights = np.zeros((self.n, self.n), dtype=np.float)
        for n1 in range(self.n):
            for n2 in range(self.n):
                edge_weights[n1, n2] = self.find_edge_weight(n1, n2, alpha) #TODO: minimizer iterations
        return edge_weights


    def solve_multihop_routing_naive(self):
        iteration_count = 0
        print("iteration: ", iteration_count, " a/t  ", self.cur_alpha, "/", self.cur_time)
        self.debug_pretty_print_init_traffic()
        while self.cur_time < self.W:

            iteration_count += 1
            remaining_time = self.W - self.cur_time
            m, alpha = self.find_M_alpha(max_duration=remaining_time)
            if m is None: break
            self.cur_time += int( alpha )
            self.cur_matching, self.cur_alpha = m, alpha
            self.matching_history.append( (self.cur_matching, self.cur_alpha) )
            self.virtual_route_flows()

            print("debug: current matching")
            print(self.cur_matching)

            print("iteration: ",iteration_count, "a/t ",self.cur_alpha,"/",self.cur_time)
            self.debug_pretty_print_current_traffic()


        print("debug: #-matching :", self.stat_matching_count,\
                " avg. time:", np.round(self.stat_matching_time/self.stat_matching_count, 3), "seconds")
        return

    def solve_multihop_routing(self):
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
                src, dest, nh, f,_,_ = val
                print("\t(", src,",", dest,")->",nh," f:",f)
        print("=======================================")
        return
    #-------end of class definition--------------#

if __name__ == '__main__':
    '''
    module test
    '''