#%%
# import Solver
from copy import deepcopy as cp

import sys
import os

# Import Admk solver for graphs
print(os.getcwd())
sys.path.append('../src/')
# sys.path.append('./src/') # running from console
from admk import Graph
from admk import MinNorm
from admk import TdensPotentialVelocity
from admk import AdmkControls
from admk import AdmkSolver

import numpy as np
from scipy.linalg import norm
#from scipy import linalg
import time as cputiming
import os
import random
import networkx as nx
import pickle
import pandas as pd
import matplotlib.cm as cmx

#%%
def get_topology(filename, nodes=None):
    df = pd.read_csv(filename, sep=',')
    if nodes is not None:
        df = df[(df['source'].isin(nodes)) & (df['target'].isin(nodes))]
    topol = np.array([df['source'].values, df['target'].values])
    return topol, df['length'].values

def get_transfer_nodes(filename):
    """
    Get the list of nodes where the transfer is non-zero
    """
    return pd.read_csv(filename, sep=',').iloc[:,0].values


def get_tranfert(filename,station_id):
    """
    Get the transfer from a station to all the nodes
    """
    df = pd.read_csv(filename, sep=',', index_col=0)
    return df.loc[str(station_id)].fillna(0).values

def node2index(nodes_list):
    """
    Given a list of nodes, return
    node2index: node2index[node_id] = i and -1 if node_id is not in nodes_list
    """
    print(nodes_list)
    n_nodes = np.amax(nodes_list)
    print(np.amin(nodes_list),n_nodes)
    node2index = - np.ones(n_nodes+1, dtype=int)*2*n_nodes
    for i, node in enumerate(nodes_list):
        node2index[node] = i
    return node2index

def run_main_multicommodity(od_matrix, origin_nodes, verbose=0):
    matrix_nodes = pd.read_csv(od_matrix_fp,
                               index_col=0).index.values
    # read topology
    topol, weight = get_topology('NetworkData/edge_attributes.csv',
                                 nodes=matrix_nodes)

    # get nodes list
    nodes_list = np.unique(topol)

    # build a map from node to index in nodes_list
    node2ind = node2index(nodes_list)

    # convert topology to be index based
    topol_index = np.array([node2ind[topol[0]], node2ind[topol[1]]])

    # Init. graph problem
    graph=Graph(topol_index)
    print('graph n_edges',graph.n_edges,'n_nodes',graph.n_nodes)

    # get the list of nodes involved in the transfer
    # transfer_nodes = get_transfer_nodes(od_matrix_fp)
    transfer_nodes = od_matrix.columns.values

    # create rhs term from a given station, say, the first one
    forcings = []
    for root_node in origin_nodes:
        # rhs = get_tranfert('NetworkData/rods_station_am_matrix_nodes.csv',
        #                    root_node)
        rhs = od_matrix.loc[root_node].fillna(0).values
        forcing = np.zeros(len(nodes_list))
        # balance the mass
        forcing[node2ind[transfer_nodes]] = rhs
        mass = forcing.sum()
        forcing[node2ind[root_node]] = -mass
        forcings.append(forcing)

    forcing = np.concatenate(forcings)

    # Init. signed incidence matrix
    incidence_matrix = graph.signed_incidence_matrix()
    incidence_matrix_transpose = incidence_matrix.transpose()

    # Init problem (same graph)
    print(incidence_matrix_transpose.size)
    problem = MinNorm(incidence_matrix_transpose, weight)

    # set problem inputs (forcing loads, powers, etc) and check
    problem = problem.set_inputs(forcing, 1.0)
    consistency = problem.check_inputs()

    # Init container for transport problem solution with
    # solution.tdens=edge conductivity
    # solution.pot=potential
    # solution.flux=conductivity * potential gradient
    # Init solver
    admk = AdmkSolver(problem)
    solution = TdensPotentialVelocity(admk.n_tdens, admk.n_pot*admk.problem.n_rhs)

    # Init solver controls
    ctrl = AdmkControls()

    # mehtod and max_iter
    ctrl.time_discretization_method = 'explicit_tdens'
    ctrl.max_iter = 200

    # deltat controls
    ctrl.deltat_control = 'expanding'
    ctrl.deltat = 1e-1
    ctrl.min_deltat = 1e-2
    ctrl.max_deltat = 5e-1

    # verbosity
    ctrl.verbose = verbose

    # solve
    ierr = admk.solve(problem, solution, ctrl)
    print('ierr=',ierr,admk.ierr_dictionary(ierr))

    # reshaping to rows are origin nodes
    solution.pot = np.reshape(solution.pot,
                              [len(origin_nodes), len(transfer_nodes)])
    # check if convergence is achieved
    return solution, problem, topol, topol_index

#%%

#
# #%%
# # Save solution as pickle
# admk_dict = {'topol': topol, 'topol_ind': topol_index, 'solution': solution,
#              'problem': problem, 'nodes': np.unique(topol)}
# pickle.dump(admk_dict, open('./results/london_mc_total.p', 'wb'))

#%%
if __name__ == "__main__":
    od_matrix_fp = 'NetworkData/rods_station_am_matrix_nodes.csv'

    od_matrix = pd.read_csv(od_matrix_fp, sep=',', index_col=0)
    od_matrix.columns = od_matrix.columns.astype(int)
    od_matrix = od_matrix.drop(columns=[153])

    od_matrix = od_matrix.loc[:, od_matrix.index]

    run_main_multicommodity(od_matrix=od_matrix,
        origin_nodes=od_matrix.columns.values[:2], verbose=0)
