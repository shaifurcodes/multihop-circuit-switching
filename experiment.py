#!/usr/bin/env python
'''
Module: matching
List of entities:
Use:

'''
from __future__ import  print_function
from multihop_matching import Multihop_Matching,  ALGO_TYPE
from traffic_generator import  Traffic_Generator

import   numpy as np

def experiment_runner():
    generate_traffic = True
    generate_route_only = True
    base_file_name = './data/synthetic/synthetic_1'

    topo_file = base_file_name+'.topology.txt'
    traffic_file = base_file_name+'.traffic.txt'
    routing_file = base_file_name+'.routing.txt'

    number_of_nodes = 10
    is_complete_graph = True #--snot compatible for sparse graph yet!!
    edge_sparsity = 100.0 #--useless unless above is False
    max_long_flow = 200
    min_long_flow = 10
    sparsity = 50. #---of traffic matrix---
    skewness = 5. #---ratio between small to large flowss
    max_hop = 1 #--as in diameter---#

    W = 300   #--window size---#
    delta = 1 #--switching delay

    algo_type = [ ALGO_TYPE.NAIVE ]

    #-------first generate synthetic data---#
    if generate_traffic:
        tg = Traffic_Generator()
        tg.generate_synthetic_traffic( base_file_name = base_file_name,
                                number_of_nodes = number_of_nodes,
                                is_complete_graph=is_complete_graph,
                                edge_sparsity=edge_sparsity,
                                max_long_flow=max_long_flow,
                                min_long_flow=min_long_flow,
                                sparsity=sparsity,
                                skewness=skewness,
                                max_hop=max_hop,
                                generate_route_only=generate_route_only)

    #----now run experiments----#
    mmatching = Multihop_Matching(W=W, delta=delta, topo_file=topo_file, traffic_file=traffic_file, routing_file=routing_file)
    for a_type in algo_type:
        mmatching.solve_multihop_routing(algo_type=a_type)

        demand_met = mmatching.find_demand_met()

        print("Demand Met: ", 100.*demand_met,"%")
        return

if __name__ == '__main__':
    np.random.seed(0)
    experiment_runner()