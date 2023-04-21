#%%
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pickle
import networkx as nx
import geopandas as gpd
import mplcyberpunk as mplcp
import matplotlib.cm as cmx

print(os.getcwd())
sys.path.append('../src/')

from admk import *
#%%
def remove_axes(ax):

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.tick_params(left = False, right = False , labelleft = False ,
                    labelbottom = False, bottom = False)

    return ax
#%%
# count_csv = pd.read_csv('./NetworkData/rods_station_am_matrix_nodes.csv')
# count_csv = count_csv.groupby(
#     ['OrigName', 'DestName'], as_index=False).agg('sum')

ntwrk_graph = pickle.load(open('../../urban-growth/LondonNetwork/core_data'
    '/network/network2011.p', 'rb'))
pos_attr = nx.get_node_attributes(ntwrk_graph, 'pos')
pos_dict = {}
for n in pos_attr.keys():
    pos_dict[n] = np.array(pos_attr[n]) * 1000
pos_dict[660] = [509348, 200000]

xvyv = pickle.load(
    open("../../urban-growth/LondonNetwork/growth_results/xvyv.p", "rb"))
rho_2011 = pickle.load(
    open("../../urban-growth/LondonNetwork/growth_results/rho_yrs.p", 'rb'))[
        2011]
matrix_nodes = pd.read_csv('NetworkData/rods_station_total_matrix_nodes.csv', 
                           index_col=0).index.values

#%%
admk_results = pickle.load(open('./results/london_test_total.p', 'rb'))
orig_topol = admk_results['topol']
topol_ind = admk_results['topol_ind']
solution = admk_results['solution']
problem = admk_results['problem']

admk_ntwrk = ntwrk_graph.subgraph(admk_results['nodes'])

bounds = np.array(list(nx.get_node_attributes(admk_ntwrk, 'pos').values())) *1000
ylim = bounds[:, 1].min() - 1000, bounds[:, 1].max() + 1000
xlim = bounds[:, 0].min() - 1000, bounds[:, 0].max() + 1000

#%%
od_edges = [tuple(e) for e in orig_topol.T[:, :]]
nx.set_node_attributes(admk_ntwrk, 
    {n: solution.pot[i]  for i, n in enumerate(matrix_nodes)}, name='potential')

nx.set_edge_attributes(admk_ntwrk,
    {tuple(e): solution.tdens[i] for i, e in enumerate(od_edges)}, 
    name='tdens') 

nx.set_edge_attributes(admk_ntwrk,
    {tuple(e): solution.flux[i] for i, e in enumerate(od_edges)}, 
    name='flux') 
#%%

plt.style.use("default")
fig, ax = plt.subplots(figsize=(6.4,3.6))
# fig = plt.figure(figsize=(10,5))

# ax = fig.add_axes([0.2, 0, 0.8, 1])
# ax.grid(b=None)
ax.patch.set_facecolor('grey')

# ax = remove_axes(ax=ax)

# ax.set_xlim(495000.0, 568000.0)
# ax.set_ylim(159000.0, 208000.0)

# ax.set_xlim(xvyv['xv'].min()*1000+5000, xvyv['xv'].max()*1000-5000)
# ax.set_ylim(xvyv['yv'].min()*1000+11000,  xvyv['yv'].max()*1000-5000)

# low_cond_edge = []
# zero_cond_edge = []
# for e in orig_topol.T[:, :]:
#     u = e[0]
#     v = e[1]
#     x = [pos_dict[u][0],pos_dict[v][0]]
#     y = [pos_dict[u][1],pos_dict[v][1]]
#     lw = np.log10(admk_ntwrk.edges[(u,v)]['tdens'])
#     if lw < 0:
#         lw=0.1
#         low_cond_edge.append(tuple(e))

#     l = Line2D(x,y,linewidth=lw, solid_capstyle='round', color='limegreen')
#     ax.add_line(l)

cond_colours = \
    np.log10(list(nx.get_edge_attributes(admk_ntwrk, 'tdens').values()))
cond = nx.draw_networkx_edges(admk_ntwrk, width=2, edge_color=cond_colours, 
    nodesize=0, pos=pos_dict, ax=ax)

# nx.draw_networkx_edges(admk_ntwrk, edgelist=low_cond_edge, width=0.5, 
#     edge_color='m', nodesize=0, pos=pos_dict, ax=ax)


sm = plt.cm.ScalarMappable(cmap='viridis', 
                           norm=plt.Normalize(vmin=np.nanmin(cond_colours), 
                                              vmax=np.nanmax(cond_colours)))
sm.set_array([])

cbar = plt.colorbar(sm, label='log$_{10}$$\mu$')
ax.set_title('Edge conductivity (log)')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_aspect('equal')
plt.tight_layout()
plt.savefig('./results/total_logcond.png', dpi=300)


#%%

plt.style.use("default")
fig, ax = plt.subplots(figsize=(6.4,3.6))

ax.patch.set_facecolor('grey')

cond_colours = \
    list(nx.get_edge_attributes(admk_ntwrk, 'tdens').values())
cond = nx.draw_networkx_edges(admk_ntwrk, width=2, edge_color=cond_colours, 
    nodesize=0, pos=pos_dict, ax=ax)

sm = plt.cm.ScalarMappable(cmap='viridis', 
                           norm=plt.Normalize(vmin=np.nanmin(cond_colours), 
                                              vmax=np.nanmax(cond_colours)))
sm.set_array([])

cbar = plt.colorbar(sm, label='log$_{10}$$\mu$')
ax.set_title('Edge conductivity')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_aspect('equal')
plt.tight_layout()
plt.savefig('./results/total_cond.png', dpi=300)
#%%
# fig = plt.figure(figsize=(10,5))

# ax = fig.add_axes([0.2, 0, 0.8, 1])
fig, ax = plt.subplots(figsize=(6.4,3.6))
# ax.grid(b=None)
ax.patch.set_facecolor('grey')

pot_colour = np.array(
    [admk_ntwrk.nodes[n]['potential'] if n in matrix_nodes else None
     for n in admk_ntwrk.nodes])
edge = nx.draw_networkx_edges(admk_ntwrk, pos=pos_dict, 
    node_size=0, ax=ax, color='k')

pot = nx.draw_networkx_nodes(admk_ntwrk, pos=pos_dict, #nodelist=matrix_nodes, 
    node_color=pot_colour, node_size=5, ax=ax, cmap='YlOrRd', )

sm2 = plt.cm.ScalarMappable(cmap='YlOrRd', 
                           norm=plt.Normalize(vmin=np.nanmin(pot_colour), 
                                              vmax=np.nanmax(pot_colour)))
sm2.set_array([])
cbar2 = plt.colorbar(sm2, label='$p_{v}$')
ax.set_title('Node potential')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_aspect('equal')
plt.tight_layout()
plt.savefig('./results/total_pot.png', dpi=300)
#%%
fig, ax = plt.subplots(figsize=(6.4,3.6))
# ax = fig.add_axes([0.2, 0, 0.8, 1])
# ax.grid(b=None)
ax.patch.set_facecolor('grey')
# ax = remove_axes(ax=ax)
ax.axis('equal')
# ax.set_xlim(495000.0, 568000.0)
# ax.set_ylim(159000.0, 208000.0)
ax.set_xlim(xlim)
ax.set_ylim(ylim)

# flux_colour = np.log10(np.array([admk_ntwrk.edges[e]['flux'] if e in od_edges else np.nan
#      for e in admk_ntwrk.edges]))
flux_colour = np.abs(list(nx.get_edge_attributes(admk_ntwrk, 'flux').values()))

nx.draw_networkx_edges(admk_ntwrk, width=2, edge_color=flux_colour,
   nodesize=0, pos=pos_dict, ax=ax)

sm = plt.cm.ScalarMappable(cmap='viridis', 
                           norm=plt.Normalize(vmin=np.nanmin(flux_colour), 
                                              vmax=np.nanmax(flux_colour)))
sm.set_array([])
cbar = plt.colorbar(sm, label='|$F_{e}$|')
ax.set_title('Passenger flux')

ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_aspect('equal')
plt.tight_layout()
# plt.savefig('./results/total_flux.png', dpi=300)