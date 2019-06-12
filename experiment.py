#!/usr/bin/env python
'''
Module: matching
List of entities:
Use:

'''
from __future__ import  print_function
from multihop_matching import Multihop_Matching
from traffic_generator import  Traffic_Generator
import cProfile
from timeit import default_timer as mytimer
import pstats

import  numpy as np

def experiment_runner():

    base_file_name = './data/synthetic/multipath_1'

    topo_file = base_file_name+'.topology.txt'
    traffic_file = base_file_name+'.traffic.txt'
    routing_file = base_file_name+'.routing.txt'

    number_of_nodes = 20
    cl, cs, nl, ns = 70, 30, 2, 6

    W = 1000   #--window size---#
    delta = 20 #--switching delay
    max_hop = 3
    multipath_factor = 1

    all_algo_types = [ 'naive',  'multipath', 'multipath_wbt', 'uppberbound', 'benchmark_1', 'benchmark_2']
    algo_type = all_algo_types[0]
    #-------first generate synthetic data---#
    tg = Traffic_Generator()
    tg.generate_synthetic_traffic(
                               base_file_name,
                               number_of_nodes = number_of_nodes,
                               W = W,
                               cl = cl,
                               cs =cs,
                               nl =nl,
                               ns =ns,
                               max_hop=max_hop,
                               multipath_factor=multipath_factor)




    #----now run experiments----#
    mmatching = Multihop_Matching(W=W,
                                  delta=delta,
                                  topo_file=topo_file,
                                  traffic_file=traffic_file,
                                  algo_type=algo_type,
                                  routing_file=routing_file)

    mmatching.solve_multihop_routing()

    demand_met = mmatching.find_demand_met()
    print("Input Parameters:\n==========================")
    print("# Nodes:", number_of_nodes)
    print("W, delta:", W, ",", delta)
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