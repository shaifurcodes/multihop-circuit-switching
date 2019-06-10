#!/usr/bin/env python
'''
Module: matching
List of entities:
Use:

'''
from __future__ import  print_function
from multihop_matching import Multihop_Matching,  ALGO_TYPE
from traffic_generator import  Traffic_Generator
import cProfile
from timeit import default_timer as mytimer
import pstats

import  numpy as np

def experiment_runner():
    generate_traffic = True
    generate_route_only = False
    base_file_name = './data/synthetic/synthetic_1'

    topo_file = base_file_name+'.topology.txt'
    traffic_file = base_file_name+'.traffic.txt'
    routing_file = base_file_name+'.routing.txt'

    number_of_nodes = 200
    is_complete_graph = True #--snot compatible for sparse graph yet!!
    edge_sparsity = 100.0 #--useless unless above is False

    cl, cs, nl, ns = 7000, 3000, 4, 12

    max_long_flow = 100
    min_long_flow = 20
    sparsity = 10. #---of traffic matrix---
    skewness = 5 #---ratio between small to large flowss
    max_intermediate_nodes = 2 #--as in diameter---#

    W = 10000   #--window size---#
    delta = 20 #--switching delay

    algo_type = [ ALGO_TYPE.NAIVE ]

    #-------first generate synthetic data---#
    if generate_traffic:
        tg = Traffic_Generator()
        tg.generate_synthetic_traffic( base_file_name = base_file_name,
                                number_of_nodes = number_of_nodes,
                                cl=cl, cs=cs, nl=nl, ns=ns,
                                is_complete_graph=is_complete_graph,
                                edge_sparsity=edge_sparsity,
                                max_long_flow=max_long_flow,
                                min_long_flow=min_long_flow,
                                sparsity=sparsity,
                                skewness=skewness,
                                max_hop=max_intermediate_nodes+1,
                                generate_route_only=generate_route_only)

    #----now run experiments----#
    mmatching = Multihop_Matching(W=W, delta=delta, topo_file=topo_file, traffic_file=traffic_file, routing_file=routing_file)
    for a_type in algo_type:
        mmatching.solve_multihop_routing(algo_type=a_type)

        demand_met = mmatching.find_demand_met()
        print("Input Parameters:\n==========================")
        print("# Nodes:", number_of_nodes)
        print("W, delta:", W, ",", delta)
        print("Long Flow Size: (",max_long_flow,",", min_long_flow,")")
        print("Traffic Sparsity: ", sparsity, "%" )
        print("Short-to-long Flow Ratio (Skewness) :",skewness)
        print("Max. Allowed Hop (Diameter):", max_intermediate_nodes+1," debug new_variable_check self.max_hop: ", mmatching.max_hop)
        print("Results:\n==========================")
        print("Demand Met: ", 100.*demand_met,"%")
        return

def parse_profile(profile_filename):
    p = pstats.Stats(profile_filename)
    p.sort_stats('calls', 'tottime').print_stats()

if __name__ == '__main__':
    np.random.seed(0)
    start_t = mytimer()

    profileCode = True
    if profileCode:
        cProfile.run('experiment_runner()', 'expProfile.cprof')
    else:
        experiment_runner()
    print('Execution time:', np.round((mytimer() - start_t), 3), "seconds")