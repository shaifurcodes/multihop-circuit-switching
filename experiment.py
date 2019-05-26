from multihop_matching import Multihop_Matching,  ALGO_TYPE

def experiment_runner():
    topo_file = './data/synthetic/topology.txt'
    traffic_file = './data/synthetic/traffic.txt'
    routing_file = './data/synthetic/routing.txt'

    W = 1000
    delta = 1

    algo_type = ALGO_TYPE.NAIVE

    mmatching = Multihop_Matching(W=W, delta=delta, topo_file=topo_file, traffic_file=traffic_file, routing_file=routing_file)
    mmatching.solve_multihop_routing(algo_type=algo_type)
    return

if __name__ == '__main__':
    experiment_runner()