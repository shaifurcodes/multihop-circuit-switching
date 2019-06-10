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
                                   cl,
                                   cs,
                                   nl,
                                   ns,
                                   is_complete_graph = True,
                                   edge_sparsity = 100.0,
                                   max_long_flow = 200,
                                   min_long_flow = 10,
                                   sparsity = 50.,
                                   skewness = 5,
                                   max_hop = 1,
                                   generate_route_only = False ):
        '''

        :param base_file_name:
        :param number_of_nodes:
        :param is_complete_graph:
        :param sparsity:
        :return:
        '''
        #----------write topology file--------
        if not generate_route_only:
            if is_complete_graph:
                with open(base_file_name+".topology.txt","w") as f:
                    f.write(str(number_of_nodes)+" C")
            else:
                pass #TODO: handle sparse graph-generation

        #-------------write traffic file-------
        if not generate_route_only:
            traffic = self.generate_sigmetric_traffic(number_of_nodes, cl, cs, nl, ns)
            np.savetxt(fname=base_file_name+".traffic.txt", X=traffic, fmt='%ld', delimiter=',')
            # traffic = np.zeros((number_of_nodes, number_of_nodes), dtype=np.longlong)
            # row, col = np.where(~np.eye(traffic.shape[0],dtype=bool))
            # n_nondiag = row.shape[0]
            # n_nonzero = int( np.round( sparsity*(number_of_nodes*number_of_nodes - number_of_nodes)/100.) )
            # nonzero_rc_indx = np.random.choice(n_nondiag, n_nonzero, replace=False)
            #
            #
            # for rc_indx in nonzero_rc_indx:
            #     i, j = row[rc_indx], col[rc_indx]
            #     traffic[i, j] = np.random.randint(1, min_long_flow, size = 1, dtype=np.longlong)
            #
            # n_long_flow = int( np.round(1.*n_nonzero*(1/(1+skewness)) ) )
            # long_rc_indx = np.random.choice(n_nonzero, n_long_flow, replace=False)
            # for lrc_indx in long_rc_indx:
            #     i, j = row[ nonzero_rc_indx[lrc_indx]  ], col[ nonzero_rc_indx[lrc_indx] ]
            #     traffic[i, j] = np.random.randint(min_long_flow, max_long_flow, size = 1, dtype=np.longlong)
            #
            # np.savetxt(fname=base_file_name+".traffic.txt", X=traffic, fmt='%ld', delimiter=',')

        #-----------generate routes-----------------
        #choose a random number between 1 and max_hop inclusive, choose random nodes
        traffic = np.loadtxt(fname=base_file_name+".traffic.txt", delimiter=',', dtype=np.longlong) #---test---#
        node_list = set( range(number_of_nodes) )
        with open(base_file_name+".routing.txt","w") as f:
            if max_hop >= 2:
                for i in range(number_of_nodes):
                    for j in range(number_of_nodes):
                        if traffic[i, j] > 0:
                            n_inter = np.random.randint(0, max_hop, size= 1, dtype=np.int)
                            if n_inter == 0:
                                continue
                            selectable_nodes = list( node_list - set([i, j]) )
                            inter_nodes = []
                            if len(selectable_nodes) <= n_inter:
                                inter_nodes = selectable_nodes
                            else:
                                inter_nodes =  np.random.choice(selectable_nodes, n_inter, replace=False)
                            file_txt = str(i)+" , "+str(j)+" : "+str(inter_nodes[0])
                            if len(inter_nodes)>1:
                                for k in inter_nodes[1:]:
                                    file_txt += ", "+str(k)
                            file_txt += "\n"
                            f.write(file_txt)

        #TODO: adapt it for sparse/non-complete graphs
        return

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
    # base_file_name = "./data/synthetic/synthetic_1"
    # number_of_nodes = 4
    tg = Traffic_Generator()
    # tg.generate_synthetic_traffic(base_file_name=base_file_name, number_of_nodes = number_of_nodes)
    n=200
    #print( tg.generate_permutation_matrix(10) )
    cl, cs, nl, ns = 7000, 3000, 4, 12
    x = tg.generate_sigmetric_traffic(n=n, cl=cl, cs=cs, nl=nl, ns=ns)
