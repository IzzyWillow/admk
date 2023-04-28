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

od_matrix_fp = 'NetworkData/rods_station_am_matrix_nodes.csv'

od_matrix = pd.read_csv(od_matrix_fp, sep=',', index_col=0)
od_matrix.columns = od_matrix.columns.astype(int)
od_matrix = od_matrix.drop(columns=[153])

od_matrix = od_matrix.loc[:, od_matrix.index]

solution, problem, topol, topol_index = \
    run_main_multicommodity(od_matrix=od_matrix,
                            origin_nodes=od_matrix.columns.values[:2],
                            verbose=0)

admk_dict = {'topol': topol, 'topol_ind': topol_index, 'solution': solution,
             'problem': problem, 'nodes': np.unique(topol)}
pickle.dump(admk_dict, open('./results/london_mc_total.p', 'wb'))
