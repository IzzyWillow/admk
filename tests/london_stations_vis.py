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

def check_edge_columns(df, ntwrk):
    """
    Checks if edges are swapped in the DataFrame header compared to the 
    representation in the networkx Graph. Swaps order if incorrect.
    """
    orig_name = []
    new_name = []
    for e in list(ntwrk.edges):
        if e in list(df.columns.values):
            pass
        elif (e[1], e[0]) in list(df.columns.values):
            orig_name.append((e[1], e[0]))
            new_name.append(e)
        else:
            raise KeyError(f'Neither {e} nor {(e[1], e[0])} in column names')            
    
    df = df.rename(columns=dict(zip(orig_name, new_name)))

    return df

def check_node_columns(df, ntwrk):
    """
    Checks all nodes are present in DataFrame
    """
    missing = []

    for n in list(ntwrk.nodes()):
        i

def plot_edges(data, ntwrk, label, single_node=None, title=None, xylim=None, 
               savepath=None):

    fig, ax = plt.subplots(figsize=(6.4,3.6))

    ax.patch.set_facecolor('grey')

    clrs = data
    nx.draw_networkx_edges(ntwrk, width=2, edge_color=clrs, 
        nodesize=0, pos=pos_dict, ax=ax)
    
    if single_node is not None:
        nx.draw_networkx_nodes(ntwrk, width=2, node_size=10, pos=pos_dict, 
                               node_color='r', nodelist=[single_node], ax=ax)

    sm = plt.cm.ScalarMappable(cmap='viridis', 
                               norm=plt.Normalize(vmin=np.nanmin(clrs), 
                                                  vmax=np.nanmax(clrs)))
    sm.set_array([])

    cbar = plt.colorbar(sm, label=label)
    if title is not None:
        ax.set_title(title)
    if xylim is not None:
        ax.set_xlim(xylim[0])
        ax.set_ylim(xylim[1])
    ax.set_aspect('equal')
    plt.tight_layout()

    if savepath is not None:
        plt.savefig(savepath, dpi=300)

def plot_nodes(data, ntwrk, label, title=None, xylim=None, savepath=None):
    
    fig, ax = plt.subplots(figsize=(6.4,3.6))
    # ax.grid(b=None)
    ax.patch.set_facecolor('grey')

    pot_colour = data
    edge = nx.draw_networkx_edges(ntwrk, pos=pos_dict, 
        node_size=0, ax=ax, color='k')

    pot = nx.draw_networkx_nodes(ntwrk, pos=pos_dict,
        node_color=pot_colour, node_size=5, ax=ax, cmap='YlOrRd', )

    sm2 = plt.cm.ScalarMappable(cmap='YlOrRd', 
                            norm=plt.Normalize(vmin=np.nanmin(pot_colour), 
                                                vmax=np.nanmax(pot_colour)))
    sm2.set_array([])
    cbar2 = plt.colorbar(sm2, label=label)
    if title is not None:
        ax.set_title('Node potential')
    if xylim is not None:
        ax.set_xlim(xylim[0])
        ax.set_ylim(xylim[1])
    ax.set_aspect('equal')
    plt.tight_layout()
    if savepath is not None:
        plt.savefig(savepath, dpi=300)
#%%
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
admk_results = pickle.load(open('./results/london_total_df.p', 'rb'))

orig_topol = admk_results['topol']
admk_ntwrk = ntwrk_graph.subgraph(admk_results['nodes'])

flux_df = check_edge_columns(df=admk_results['flux_df'], ntwrk=admk_ntwrk)
pot_df = admk_results['pot_df']
conduct_df = check_edge_columns(df=admk_results['conduct_df'], ntwrk=admk_ntwrk)

bounds = np.array(list(nx.get_node_attributes(admk_ntwrk, 'pos').values())) *1000
ylim = bounds[:, 1].min() - 1000, bounds[:, 1].max() + 1000
xlim = bounds[:, 0].min() - 1000, bounds[:, 0].max() + 1000
#%%
od_edges = [tuple(e) for e in orig_topol.T[:, :]]

orig_node = 52

log_cond_vals = np.log10(
    conduct_df.loc[orig_node, list(admk_ntwrk.edges)].values)
plot_edges(data=log_cond_vals, ntwrk=admk_ntwrk, label='log$_{10}$$\mu$', 
           title='Edge conductivity (log)', savepath=None, xylim=(xlim, ylim),
           single_node=orig_node)

cond_vals = conduct_df.loc[orig_node, list(admk_ntwrk.edges)].values
plot_edges(data=cond_vals, ntwrk=admk_ntwrk, label='$\mu$', xylim=(xlim, ylim),
           title='Edge conductivity', savepath=None, single_node=orig_node)

pot_vals = pot_df.loc[orig_node, list(admk_ntwrk.nodes)]
plot_nodes(data=pot_vals, ntwrk=admk_ntwrk, label='$p_{v}$', xylim=(xlim, ylim),
           title='Node potential')#, single_node=orig_node)

flux_vals = flux_df.loc[orig_node, list(admk_ntwrk.edges)].values
plot_edges(data=flux_vals, ntwrk=admk_ntwrk, label='|$F_{e}$|', 
           title='Passenger flux', savepath=None, xylim=(xlim, ylim), 
           single_node=orig_node)

#%%
od_edges = [tuple(e) for e in orig_topol.T[:, :]]

# orig_node = 52
# log_cond_vals = np.log10(
#     conduct_df.loc[orig_node, list(admk_ntwrk.edges)].values)


plot_edges(data=np.log10(conduct_df[list(admk_ntwrk.edges)].sum(axis=0).values), 
           ntwrk=admk_ntwrk, label='log$_{10}$$\Sigma_{e}\mu$', 
           title='Edge conductivity (log)', savepath=None, xylim=(xlim, ylim))

# cond_vals = conduct_df.loc[orig_node, list(admk_ntwrk.edges)].values
plot_edges(data=conduct_df[list(admk_ntwrk.edges)].sum(axis=0).values, 
           ntwrk=admk_ntwrk, label='$\Sigma_{e}\mu$', xylim=(xlim, ylim),
           title='Edge conductivity', savepath=None)

# # pot_vals = pot_df.loc[orig_node, list(admk_ntwrk.nodes)]
# plot_nodes(data=pot_df.sum(axis=0), ntwrk=admk_ntwrk, label='$p_{v}$', xylim=(xlim, ylim),
#            title='Node potential')#, single_node=orig_node)

# flux_vals = flux_df.loc[orig_node, list(admk_ntwrk.edges)].values
plot_edges(data=flux_df[list(admk_ntwrk.edges)].sum(axis=0).values, 
           ntwrk=admk_ntwrk, label='$\Sigma_{e}|F_{e}$|', xylim=(xlim, ylim),
           title='Passenger flux', savepath=None)
