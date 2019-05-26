#!/usr/bin/env python
'''
Module: matching
List of entities:
Use:

'''
from __future__ import  print_function
from multihop_matching import Multihop_Matching,  ALGO_TYPE

def experiment_runner():
    topo_file = './data/synthetic/topology.txt'
    traffic_file = './data/synthetic/traffic.txt'
    routing_file = './data/synthetic/routing.txt'

    W = 40
    delta = 1

    algo_type = ALGO_TYPE.NAIVE

    mmatching = Multihop_Matching(W=W, delta=delta, topo_file=topo_file, traffic_file=traffic_file, routing_file=routing_file)
    mmatching.solve_multihop_routing(algo_type=algo_type)
    demand_met = mmatching.find_demand_met()
    print("Demand Met: ",100.*demand_met,"%")
    return

if __name__ == '__main__':
    experiment_runner()