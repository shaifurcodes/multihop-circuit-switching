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
import random

def experiment_runner():

    base_file_name = './data/synthetic/multipath_1'

    number_of_nodes = 64
    cl, cs, nl, ns = 70, 30, 2, 6

    W = 10000
    delta = 20
    max_hop = 3
    multipath_factor = 3

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
    mmatching = Multihop_Matching( W=W,
                                  delta=delta,
                                  algo_type=algo_type,
                                  base_file_path=base_file_name )

    mmatching.run_multihop_matching_experiment()

    demand_met = mmatching.find_demand_met()
    print("Input Parameters:\n==========================")
    print("Algo Type:", mmatching.cur_algo_type)
    print("# Nodes:", number_of_nodes)
    print("W, delta:", W, ",", delta)
    print("Max Hop:", mmatching.max_hop, "  Max Path: ", mmatching.multipath_factor)
    print("Results:\n==========================")
    print("Demand Met: ", 100.*demand_met,"%")
    print("Results:\n==========================")

    print("Debug: overflow amount:", mmatching.debug_overflow_amount)
    return

def parse_profile(profile_filename):
    p = pstats.Stats(profile_filename)
    p.sort_stats('calls', 'tottime').print_stats()

if __name__ == '__main__':
    np.random.seed(11739)
    random.seed(11739)

    start_t = mytimer()

    profileCode = False
    if profileCode:
        cProfile.run('experiment_runner()', 'expProfile.cprof')
    else:
        experiment_runner()
    print('Execution time:', np.round((mytimer() - start_t), 3), "seconds")