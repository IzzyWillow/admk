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

import pandas as pd

#%%
def get_topology(filename):
    df = pd.read_csv(filename, sep=',')
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
    df = pd.read_csv(filename, sep=',')
    return df[str(station_id)].fillna(0).values

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

def test_main(verbose=0):
    # read topology
    topol, weight = get_topology('NetworkData/edge_attributes.csv')

    # get nodes list
    nodes_list = np.unique(topol)
    
    # build a map from node to index in nodes_list
    node2ind = node2index(nodes_list)
    
    # convert topology to be index based
    topol_index = np.array([node2ind[topol[0]], node2ind[topol[1]]])
    
    # get the list of nodes involved in the transfer
    transfer_nodes = get_transfer_nodes('NetworkData/rods_station_am_matrix_nodes.csv')

    # create rhs term from a given station, say, the first one
    rhs = get_tranfert('NetworkData/rods_station_am_matrix_nodes.csv', transfer_nodes[0])
    forcing = np.zeros(len(nodes_list))
    # balance the mass
    forcing[node2ind[transfer_nodes]] = rhs
    mass = forcing.sum()
    forcing[node2ind[transfer_nodes[0]]] = -mass

    print('f', forcing.size)
    print('w', weight.size)
    print('inc',topol_index.shape)
   

    # Init. graph problem
    graph=Graph(topol_index)
    print('graph size',graph.n_edges,graph.n_nodes)
    
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
    solution = TdensPotentialVelocity(graph.n_edges,graph.n_nodes)

    # Init solver
    admk = AdmkSolver()

    # Init solver controls
    ctrl = AdmkControls()
    
    # mehtod and max_iter
    ctrl.time_discretization_method = 'explicit_tdens'
    ctrl.max_iter = 1000
    
    # deltat controls
    ctrl.deltat_control = 'fixed'
    ctrl.deltat = 1e-1
    ctrl.min_deltat = 1e-1
    ctrl.max_deltat = 1e4
    
    # verbosity
    ctrl.verbose = verbose
    
    # solve
    ierr = admk.solve(problem, solution, ctrl)

    print('ierr = ' + str(ierr))
    # check if convergence is achieved
    return 0

def to_networkx(arr, tdens, pot, flux=None):
    # arr: array of edges 
    # tdens: potential velocity of edges

    # for testing only -----------------------
    for r in arr:
        numbers = list(range(1,7))
        numbers.remove(r[0])
        r[1] = random.choice(numbers)
        numbers.remove(r[1])
        r[2] = random.choice(numbers)
    # ----------------------------------------

    G = nx.Graph()

    G.add_nodes_from(np.unique(arr))

    edges = []
    for r in arr:
        for i in r[1:]:
            edges.append((r[0], i))

    G.add_edges_from(edges)

    nx.set_node_attributes(G, pot, name='potential')

    nx.set_edge_attributes(
        G, {e: tdens[i] for i, e in enumerate(edges)}, name='tdens')    
    
    if flux is not None:
        nx.set_edge_attributes(
            G, {e: flux[i] for i, e in enumerate(edges)}, name='flux')    
    
    for e in edges:
    
#%%
# Taking the functionality out of test_main() to test it
# read topology
topol, weight = get_topology('NetworkData/edge_attributes.csv')

# get nodes list
nodes_list = np.unique(topol)

# build a map from node to index in nodes_list
node2ind = node2index(nodes_list)

# convert topology to be index based
topol_index = np.array([node2ind[topol[0]], node2ind[topol[1]]])

# get the list of nodes involved in the transfer
transfer_nodes = get_transfer_nodes(
    'NetworkData/rods_station_am_matrix_nodes.csv')

# create rhs term from a given station, say, the first one
rhs = get_tranfert('NetworkData/rods_station_am_matrix_nodes.csv', transfer_nodes[0])

forcing = np.zeros(len(nodes_list))
# balance the mass
forcing[node2ind[transfer_nodes]] = rhs
mass = forcing.sum()  # total passengers in flow
forcing[node2ind[transfer_nodes[0]]] = -mass

print('f', forcing.size)
print('w', weight.size)
print('inc',topol_index.shape)

#%%
# Init. graph problem
graph=Graph(topol_index)
print('graph size',graph.n_edges,graph.n_nodes)

# Init. signed incidence matrix
# Entry is +1 if node is starting point, -1 if end point
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
solution = TdensPotentialVelocity(graph.n_edges,graph.n_nodes)

# Init solver
admk = AdmkSolver()

# Init solver controls
ctrl = AdmkControls()

# mehtod and max_iter
ctrl.time_discretization_method = 'explicit_tdens'
ctrl.max_iter = 1000

# deltat controls
ctrl.deltat_control = 'fixed'
ctrl.deltat = 1e-1
ctrl.min_deltat = 1e-1
ctrl.max_deltat = 1e4

# verbosity
ctrl.verbose = 0
#%%
# solve
ierr = admk.solve(problem, solution, ctrl)

# # check if convergence is achieved
# return 0
#%%
if __name__ == "__main__":
    sys.exit(test_main(2))
