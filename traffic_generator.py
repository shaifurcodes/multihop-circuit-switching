#!/usr/bin/env python
'''
Module: matching
List of entities:
Use:

'''
from __future__ import  print_function

__author__ = "Md Shaifur Rahman"
__email__ = "mdsrahman@cs.stonybrook.edu"
__date__ = "May 26, 2019"
__license__ = "GPL"
__copyright__ = "Copyright 2019, WINGS Lab, Stony Brook University"


#---start here--->
import numpy as np
class Traffic_Generator(object):
    def __init__(self):
        '''

        '''
        return

    def generate_synthetic_traffic(self, base_file_name,
                                   number_of_nodes,
                                   W,
                                   cl,
                                   cs,
                                   nl,
                                   ns,
                                   max_hop = 1,
                                   multipath_factor = 1):
        '''

        :param base_file_name:
        :param number_of_nodes:
        :param W:
        :param cl:
        :param cs:
        :param nl:
        :param ns:
        :param max_hop:
        :param multipath_factor:
        :return:
        '''

        #----------write topology file--------
        with open(base_file_name+".topology.txt","w") as f:
            f.write(str(number_of_nodes)+" C")


        #-------------write traffic file-------
        cl_val = int( round(W*cl/100.0,2))
        cs_val = int( round( W*cs/100.,2))
        nl_val = int( round(number_of_nodes*nl/100.0))
        ns_val = int( round(number_of_nodes*ns/100.0))
        traffic = self.generate_sigmetric_traffic(number_of_nodes, cl_val, cs_val, nl_val, ns_val)
        np.savetxt(fname=base_file_name+".traffic.txt", X=traffic, fmt='%ld', delimiter=',')

        #-----------generate routes-----------------
        traffic_x, traffic_y = traffic.shape
        route =  [ [] for i in range(number_of_nodes) ]
        for i in range(traffic_x):
            route.append([])

        with open(base_file_name+".routing.txt","w") as f:
            for i in range(traffic_x):
                for j in range(traffic_y):
                    if traffic[i][j]>0:
                        path_list = self.get_paths(src=i, dest=j,
                                                   number_of_nodes=number_of_nodes,
                                                   max_hop=max_hop,
                                                   number_of_paths=multipath_factor)
                        if len(path_list) <= 0:
                            file_txt = str(i) + " , " + str(j) + " : "+'\n'
                            f.write(file_txt)
                        else:
                            for path in path_list:
                                file_txt = str(i) + " , " + str(j) + " : "
                                if len(path)>0:
                                    file_txt += str(path[0])
                                for k in path[1:]:
                                    file_txt += " , "+str(k)
                                file_txt += '\n'
                                f.write(file_txt)
        return

    def get_paths(self, src, dest, number_of_nodes, max_hop, number_of_paths ):
        '''
        :param src:
        :param dest:
        :param max_hop:
        :param number_of_paths:
        :return:
        '''
        max_retry = number_of_nodes
        path_list = []
        if max_hop <= 1:
            return [path_list]

        candidate_nodes = list(set( range(number_of_nodes) ) - set([src, dest]) )
        path_count = 0
        retry = 0
        while path_count < (number_of_paths):
            cur_path = []
            n_inter = np.random.randint(0, max_hop, size=1, dtype=np.int)
            inter_nodes = []
            if n_inter != 0:
                if len(candidate_nodes) <= n_inter:
                    inter_nodes = candidate_nodes
                else:
                    inter_nodes = list(np.random.choice(candidate_nodes, n_inter, replace=False))
            if inter_nodes in path_list:
                retry += 1
                if retry >= max_retry:
                    break
                else:
                    continue
            path_list.append(inter_nodes)
            path_count += 1
        return path_list

    def generate_permutation_matrix(self, n):
        '''

        :param n:
        :return:
        '''
        permutated_matrix =  np.zeros((n, n), dtype=np.int)
        shuffled_indices = np.arange(0, n)
        np.random.shuffle( shuffled_indices )
        for indx, v in np.ndenumerate(shuffled_indices):
            i, j = indx[0], v
            permutated_matrix[i, j] = 1

        return  permutated_matrix

    def generate_sigmetric_traffic(self, n, cl, cs, nl, ns):
        '''

        :param n:
        :return:
        '''
        T_l = np.zeros((n, n), dtype=np.float)
        T_s = np.zeros((n, n), dtype=np.float)
        T = np.zeros((n, n), dtype=np.float)

        if nl>0 and cl >0:
            for i in range(nl):
                T_l += self.generate_permutation_matrix(n)
            T_l = (cl / nl) * T_l

        if ns>0 and cs >0:
            for i in range(ns):
                T_s += self.generate_permutation_matrix(n)
            T_s = (cs / ns) * T_s

        T =  T_l+T_s
        #add noise
        for i in range(n):
            for j in range(n):
                if(T[i, j]>0):
                    if i==j:
                        T[i, j] = 0
                    else:
                        noise = np.round(np.random.normal(0, 0.3*10000/100), 0)
                        T[i, j] += noise

        return T.astype(int)


if __name__ == '__main__':
    base_file_name = "./data/synthetic/multipath_1"
    number_of_nodes = 10
    W = 10000
    cl, cs, nl, ns = 70, 30, 2, 6
    max_hop = 5
    multipath_factor = 3

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
