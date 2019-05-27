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
                                   is_complete_graph = True,
                                   edge_sparsity = 100.0,
                                   max_long_flow = 200,
                                   min_long_flow = 10,
                                   sparsity = 50.,
                                   skewness = 5,
                                   max_hop = 3,
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
            traffic = np.zeros((number_of_nodes, number_of_nodes), dtype=np.longlong)
            row, col = np.where(~np.eye(traffic.shape[0],dtype=bool))
            n_nondiag = row.shape[0]
            n_nonzero = int( np.round( sparsity*(number_of_nodes*number_of_nodes - number_of_nodes)/100.) )
            nonzero_rc_indx = np.random.choice(n_nondiag, n_nonzero, replace=False)


            for rc_indx in nonzero_rc_indx:
                i, j = row[rc_indx], col[rc_indx]
                traffic[i, j] = np.random.randint(1, min_long_flow, size = 1, dtype=np.longlong)

            n_long_flow = int( np.round(1.*n_nonzero*(1/(1+skewness)) ) )
            long_rc_indx = np.random.choice(n_nonzero, n_long_flow, replace=False)
            for lrc_indx in long_rc_indx:
                i, j = row[ nonzero_rc_indx[lrc_indx]  ], col[ nonzero_rc_indx[lrc_indx] ]
                traffic[i, j] = np.random.randint(min_long_flow, max_long_flow, size = 1, dtype=np.longlong)

            np.savetxt(fname=base_file_name+".traffic.txt", X=traffic, fmt='%ld', delimiter=',')

        #-----------generate routes-----------------
        #choose a random number between 1 and max_hop inclusive, choose random nodes
        traffic = np.loadtxt(fname=base_file_name+".traffic.txt", delimiter=',', dtype=np.longlong) #---test---#
        node_list = set( range(number_of_nodes) )
        with open(base_file_name+".routing.txt","w") as f:
            for i in range(number_of_nodes):
                for j in range(number_of_nodes):
                    if traffic[i, j] > 0:
                        n_inter = np.random.randint(0, max_hop+1, size= 1, dtype=np.int)
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

if __name__ == '__main__':
    base_file_name = "./data/synthetic/synthetic_1"
    number_of_nodes = 4
    tg = Traffic_Generator()
    tg.generate_synthetic_traffic(base_file_name=base_file_name, number_of_nodes = number_of_nodes)