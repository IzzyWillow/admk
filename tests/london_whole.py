"""
Additive estimation of entire network flux
"""

#%%
# import Solver 
from copy import deepcopy as cp
from copy import copy
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
import matplotlib.pyplot as plt
import math

#%%
class Recorder:
    def __init__(self) -> None:
        self.rec_dict = {}

    def rec_entry(self, k, v):
        self.rec_dict[k] = v

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
    #THIS SHOULD BE THE OTHER WAY ROUND FOR THE CURRENT OD MATRIX  
    df = pd.read_csv(filename, sep=',')
    return df[str(station_id)].fillna(0).values

def node2index(nodes_list):
    """
    Given a list of nodes, return
    node2index: node2index[node_id] = i and -1 if node_id is not in nodes_list
    """
    # print(nodes_list)
    n_nodes = np.amax(nodes_list)
    # print(np.amin(nodes_list),n_nodes)
    node2index = - np.ones(n_nodes+1, dtype=int)*2*n_nodes
    for i, node in enumerate(nodes_list):
        node2index[node] = i
    return node2index
    
def run_main(od_matrix, station_id, verbose=0, rec_forcing=None,
             rec_ind_mat=None):
    # read topology
    topol, weight = get_topology('NetworkData/edge_attributes.csv', 
                                  nodes=od_matrix.index.values)

    # get nodes list
    nodes_list = np.unique(topol)
    
    # build a map from node to index in nodes_list
    node2ind = node2index(nodes_list)
    
    # convert topology to be index based
    topol_index = np.array([node2ind[topol[0]], node2ind[topol[1]]])
    
    # get the list of nodes involved in the transfer
    transfer_nodes = od_matrix.columns.values

    # create rhs term from a given station, say, the first one
    rhs = od_matrix.loc[station_id].fillna(0).values
    forcing = np.zeros(len(nodes_list))
    # balance the mass: outlet = inlet
    forcing[node2ind[transfer_nodes]] = rhs
    mass = forcing.sum()
    stn_ind = np.argwhere(transfer_nodes == station_id)[0][0]
    forcing[node2ind[transfer_nodes[stn_ind]]] = -mass
    if rec_forcing is not None:
         rec_forcing.rec_entry(k=s_id, v=copy(forcing))

    # print('f', forcing.size)
    # print('w', weight.size)
    # print('inc',topol_index.shape)

    # Init. graph problem
    graph=Graph(topol_index)
    print('graph size',graph.n_edges,graph.n_nodes)
    
    # Init. signed incidence matrix
    incidence_matrix = graph.signed_incidence_matrix()
    incidence_matrix_transpose = incidence_matrix.transpose()
    
    # Init problem (same graph)
    # print(incidence_matrix_transpose.size)
    problem = MinNorm(incidence_matrix_transpose, weight)
    if rec_ind_mat is not None:
         rec_ind_mat.rec_entry(k=s_id, v=incidence_matrix)

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
    print('ierr=', ierr, admk.ierr_dictionary(ierr))
    # check if convergence is achieved
    return solution, problem

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

    # edges = []
    # for r in arr:
    #     for i in r[1:]:
    #         edges.append((r[0], i))

    G.add_edges_from(edges)

    nx.set_node_attributes(G, pot, name='potential')

    nx.set_edge_attributes(
        G, {e: tdens[i] for i, e in enumerate(edges)}, name='tdens')    
    
    if flux is not None:
        nx.set_edge_attributes(
            G, {e: flux[i] for i, e in enumerate(edges)}, name='flux')    
    
    # for e in edges:
    
#%%
od_matrix_fp = 'NetworkData/rods_station_total_matrix_nodes.csv'
od_matrix = pd.read_csv(od_matrix_fp, sep=',', index_col=0)
od_matrix.columns = od_matrix.columns.astype(int)

od_matrix = od_matrix.drop(columns=[153])

stations = od_matrix.index.values
topol, weight = get_topology('NetworkData/edge_attributes.csv', 
                             nodes=od_matrix.index.values)

flux = np.zeros(topol.shape[1])
pot = np.zeros_like(stations).astype(float)
conduct = np.zeros(topol.shape[1])
flux = {}
pot = {}
conduct = {}
#%%
forc_rec = Recorder()
ind_rec = Recorder()
start_all = cputiming.perf_counter()
for s_id in stations:
    print(s_id)
    solution, problem = run_main(od_matrix=od_matrix, station_id=s_id,)
    flux[s_id] = np.abs(solution.flux)
    pot[s_id] = solution.pot
    conduct[s_id] = np.abs(solution.tdens)
end_all = cputiming.perf_counter()

# elapsed times
e_time = end_all-start_all  # Roughly 6m10s
print(f"Elapsed time: {math.floor(e_time/60)}:{int(e_time%60):02}")

#%%
flux_df = pd.DataFrame.from_dict(flux, orient='index', 
                                 columns=list(zip(topol[0], topol[1])))
pot_df = pd.DataFrame.from_dict(pot, orient='index', 
                                columns=od_matrix.columns.values)
conduct_df = pd.DataFrame.from_dict(conduct, orient='index', 
                                    columns=list(zip(topol[0], topol[1])))
#%%
# Save solution as pickle
admk_dict = {'topol': topol, 'flux_df': flux_df, 'pot_df': pot_df, 
             'conduct_df': conduct_df, 'nodes': np.unique(topol)}

pickle.dump(admk_dict, open('./results/london_total_df.p', 'wb'))

#%%
if __name__ == "__main__":
    sys.exit(run_main(2))
