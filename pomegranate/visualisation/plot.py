# Custom plot functions 


import logging
from itertools import count
from typing import Dict, List, Optional, Tuple, Union

import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D

from graphein.utils.utils import protein_letters_3to1_all_caps as aa3to1

'''
Modified from graphein.protein.visualisation
'''
def motif_plot_distance_matrix(
    g: Optional[nx.Graph],
    psite: Union[int, str],
    dist_mat: Optional[np.ndarray] = None,
    use_plotly: bool = True,
    aa_order: Optional[str] = "hydro",
    reverse_order: bool = False,
    title: Optional[str] = None,
    show_residue_labels: bool = True,
) -> go.Figure:
    """Plots a distance matrix of the graph.

    :param g: NetworkX graph containing a distance matrix as a graph attribute (``g.graph['dist_mat']``).
    :type g: nx.Graph, optional
    :param psite: Location of the phosphorylation site.
    :type psite: Union[int, str]
    :param dist_mat: Distance matrix to plot. If not provided, the distance matrix is taken from the graph. Defaults to ``None``.
    :type dist_mat: np.ndarray, optional
    :param use_plotly: Whether to use ``plotly`` or ``seaborn`` for plotting. Defaults to ``True``.
    :type use_plotly: bool
    :param title: Title of the plot.Defaults to ``None``.
    :type title: str, optional
    :param aa_order: Method used to order residues on the axes.  Defaults to ``sequence`` order.
    :type aa_order: str, optional 
    :param reverse_order: Reverse the ordering of residues on the axes.  Defaults to ``False``.
    :show_residue_labels: Whether to show residue labels on the plot. Defaults to ``True``.
    :type show_residue_labels: bool
    :raises: ValueError if neither a graph ``g`` or a ``dist_mat`` are provided.
    :return: Plotly figure.
    :rtype: px.Figure
    """
    if g is None and dist_mat is None:
        raise ValueError("Must provide either a graph or a distance matrix.")

    if dist_mat is None:
        dist_mat = g.graph["distmat"] 
        
        
    # Phospho site 
    if isinstance(psite, str):
        try:
            psite = int(psite.split(':')[-1])
        except ValueError:
            raise ValueError("Specified phospho site isn't in correct format.")  
            
    # TODO Check if graph is original (i.e. size of dist_mat == size of graph)
    def get_indexes(node_list):
        indexes = []
        for n in node_list:
            indexes.append(int(n.split(":")[2]) - 1)    # indexing from 0
        return indexes  
        
        # TODO 
        # sort indexes in - distance order, - seq position order           
    	
    if g is not None:    
 
        x_range = list(g.nodes)
        
        # Sort nodes by ordering

        # TODO: store sorting function inside a dict so can be used like a case-switch thing
        if aa_order == 'seq':
            x_range = sorted(x_range, key=lambda x: int(x.split(':')[-1]))
        elif aa_order == 'euclidean':
            # sort ascending distance
            pass
        elif aa_order == 'hydro':
            # hydrophobicity ascending
            ordering = "IVLFCMAWGTSYPHNDQEKR"
            x_range = sorted(x_range, key=lambda res: [ordering.index(a) for a in aa3to1(res.split(':')[1])]) 

        else: 
            raise ValueError(f"'{aa_order}' isn't a valid axis ordering.")
         
        
        # sequence order (ascending)


        
        y_range = x_range.copy()
        
        # TODO sort by distance
        
        # add phospho label to selected site
        for i in range(len(x_range)):
            s = x_range[i]
            if psite == int(s.split(':')[-1]):
                break
        
        x_range[i] = f"{x_range[i]} [P]"
        y_range[i] = f"[P] {y_range[i]}"
        
        
        
        if not title:
            title = g.graph["name"] + " - Distance Matrix"
    else:
        x_range = list(range(dist_mat.shape[0]))
        y_range = list(range(dist_mat.shape[1]))
        if not title:
            title = "Distance matrix"

    if use_plotly:
        fig = px.imshow(
            dist_mat,
            x=x_range,
            y=y_range,
            labels=dict(color="Distance"),
            title=title,
            template="plotly_dark",
            color_continuous_scale="viridis_r",
            
        )
    else:
        if show_residue_labels:
            tick_labels = x_range
        else:
            tick_labels = []
        fig = sns.heatmap(
            dist_mat, xticklabels=tick_labels, yticklabels=tick_labels
        ).set(title=title)

    return fig
    


