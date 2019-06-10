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

    number_of_nodes = 50
    is_complete_graph = True #--snot compatible for sparse graph yet!!
    edge_sparsity = 100.0 #--useless unless above is False

    cl, cs, nl, ns = 70, 30, 2, 6

    # max_long_flow = 100
    # min_long_flow = 20
    # sparsity = 10. #---of traffic matrix---
    # skewness = 5 #---ratio between small to large flowss
    max_intermediate_nodes = 0 #--as in diameter---#

    W = 10000   #--window size---#
    delta = 20 #--switching delay

    algo_type = [ "multipath_wbt"]

    #-------first generate synthetic data---#
    if generate_traffic:
        tg = Traffic_Generator()
        tg.generate_synthetic_traffic( base_file_name = base_file_name,
                                number_of_nodes = number_of_nodes,
                                W=W,
                                cl=cl, cs=cs, nl=nl, ns=ns,
                                is_complete_graph=is_complete_graph,
                                # edge_sparsity=edge_sparsity,
                                # max_long_flow=max_long_flow,
                                # min_long_flow=min_long_flow,
                                # sparsity=sparsity,
                                # skewness=skewness,
                                max_hop=max_intermediate_nodes+1,
                                generate_route_only=generate_route_only)

    algo_type = "multipath"

    #----now run experiments----#
    mmatching = Multihop_Matching(W=W, delta=delta, topo_file=topo_file, traffic_file=traffic_file, algo_type=algo_type, routing_file=routing_file)
    mmatching.solve_multihop_routing()

    demand_met = mmatching.find_demand_met()
    print("Input Parameters:\n==========================")
    print("# Nodes:", number_of_nodes)
    print("W, delta:", W, ",", delta)
    #print("Long Flow Size: (",max_long_flow,",", min_long_flow,")")
    #print("Traffic Sparsity: ", sparsity, "%" )
    #print("Short-to-long Flow Ratio (Skewness) :",skewness)
    print("Max. Allowed Hop (Diameter):", mmatching.max_hop)
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