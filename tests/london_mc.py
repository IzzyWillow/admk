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
from london_stations import *

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

    # Init. graph problem
    graph=Graph(topol_index)
    print('graph n_edges',graph.n_edges,'n_nodes',graph.n_nodes)

    # get the list of nodes involved in the transfer
    transfer_nodes = od_matrix.columns.values

    # create rhs term from a given station, say, the first one
    forcings = []
    for root_node in transfer_nodes[[0,1]]:
        rhs = get_tranfert('NetworkData/rods_station_am_matrix_nodes.csv', root_node)
        forcing = np.zeros(len(nodes_list))
        # balance the mass
        forcing[node2ind[transfer_nodes]] = rhs
        mass = forcing.sum()
        forcing[node2ind[root_node]] = -mass
        forcings.append(forcing)
    forcing = np.concatenate(forcings)
    # print('f', forcing.size)
    # print('w', weight.size)
    # print('inc',topol_index.shape)

    # Init. signed incidence matrix
    incidence_matrix = graph.signed_incidence_matrix()
    incidence_matrix_transpose = incidence_matrix.transpose()

    # Init problem (same graph)
    # print(incidence_matrix_transpose.size)
    problem = MinNorm(incidence_matrix_transpose, weight)

    # if rec_ind_mat is not None:
    #     rec_ind_mat.rec_entry(k=s_id, v=incidence_matrix)

    # set problem inputs (forcing loads, powers, etc) and check
    problem = problem.set_inputs(forcing, 1.0)
    consistency = problem.check_inputs()

    # Init container for transport problem solution with
    # solution.tdens=edge conductivity
    # solution.pot=potential
    # solution.flux=conductivity * potential gradient
    # Init solver
    admk = AdmkSolver(problem)
    solution = TdensPotentialVelocity(admk.n_tdens,
                                      admk.n_pot*admk.problem.n_rhs)

    # Init solver controls
    ctrl = AdmkControls()

    # mehtod and max_iter
    ctrl.time_discretization_method = 'explicit_tdens'
    ctrl.max_iter = 200

    # deltat controls
    ctrl.deltat_control = 'fixed'
    ctrl.deltat = 1e-1
    ctrl.min_deltat = 1e-2
    ctrl.max_deltat = 5e-1

    # verbosity
    ctrl.verbose = verbose

    # solve
    ierr = admk.solve(problem, solution, ctrl)

    print('ierr = ' + str(ierr))
    print('ierr=', ierr, admk.ierr_dictionary(ierr))
    # check if convergence is achieved
    return solution, problem

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
