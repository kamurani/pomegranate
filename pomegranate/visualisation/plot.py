"""Custom plot functions.  Modified from graphein.protein.visualisation"""


from __future__ import annotations

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

from graphein.protein.subgraphs import extract_k_hop_subgraph
from graphein.utils.utils import import_message

# imports from visualisation.py
from graphein.protein.visualisation import colour_nodes, colour_edges


# added imports
from graphein.utils.utils import protein_letters_3to1_all_caps as aa3to1

'''
Modified from graphein.protein.visualisation by Cam M
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
    colour: Optional[str] = 'viridis_r'
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
            color_continuous_scale=colour,
            
        )
    else:
        if show_residue_labels:
            tick_labels = x_range
        else:
            tick_labels = []
        fig = sns.heatmap(
            dist_mat, xticklabels=tick_labels, yticklabels=tick_labels, cmap=colour
        ).set(title=title)

    return fig




'''
Modified from graphein.protein.visualisation
'''
def motif_plot_plotly_protein_structure_graph(
    G: nx.Graph,
    plot_title: Optional[str] = None,
    figsize: Tuple[int, int] = (620, 650),
    node_alpha: float = 0.7,
    node_size_min: float = 20.0,
    node_size_multiplier: float = 20.0,
    label_node_ids: bool = True,
    node_colour_map=plt.cm.plasma,
    edge_color_map=plt.cm.plasma,
    colour_nodes_by: str = "degree",
    colour_edges_by: Optional[str] = "kind",
) -> go.Figure:
    """
    Plots protein structure graph using plotly.

    :param G:  nx.Graph Protein Structure graph to plot
    :type G: nx.Graph
    :param plot_title: Title of plot, defaults to ``None``.
    :type plot_title: str, optional
    :param figsize: Size of figure, defaults to ``(620, 650)``.
    :type figsize: Tuple[int, int]
    :param node_alpha: Controls node transparency, defaults to ``0.7``.
    :type node_alpha: float
    :param node_size_min: Specifies node minimum size. Defaults to ``20.0``.
    :type node_size_min: float
    :param node_size_multiplier: Scales node size by a constant. Node sizes reflect degree. Defaults to ``20.0``.
    :type node_size_multiplier: float
    :param label_node_ids: bool indicating whether or not to plot ``node_id`` labels. Defaults to ``True``.
    :type label_node_ids: bool
    :param node_colour_map: colour map to use for nodes. Defaults to ``plt.cm.plasma``.
    :type node_colour_map: plt.cm
    :param edge_color_map: colour map to use for edges. Defaults to ``plt.cm.plasma``.
    :type edge_color_map: plt.cm
    :param colour_nodes_by: Specifies how to colour nodes. ``"degree"``, ``"seq_position"`` or a node feature.
    :type colour_nodes_by: str
    :param colour_edges_by: Specifies how to colour edges. Currently only ``"kind"`` or ``None`` are supported.
    :type colour_edges_by: Optional[str]
    :returns: Plotly Graph Objects plot
    :rtype: go.Figure
    """

    # Get Node Attributes
    pos = nx.get_node_attributes(G, "coords")

    # Get node colours
    node_colors = colour_nodes(
        G, colour_map=node_colour_map, colour_by=colour_nodes_by
    )
    edge_colors = colour_edges(
        G, colour_map=edge_color_map, colour_by=colour_edges_by
    )

    # 3D network plot
    x_nodes = []
    y_nodes = []
    z_nodes = []
    node_sizes = []
    node_labels = []

    # Loop on the pos dictionary to extract the x,y,z coordinates of each node
    for i, (key, value) in enumerate(pos.items()):
        x_nodes.append(value[0])
        y_nodes.append(value[1])
        z_nodes.append(value[2])
        node_sizes.append(node_size_min + node_size_multiplier * G.degree[key])

        if label_node_ids:
            node_labels.append(list(G.nodes())[i])

    nodes = go.Scatter3d(
        x=x_nodes,
        y=y_nodes,
        z=z_nodes,
        mode="markers",
        marker={
            "symbol": "circle",
            "color": node_colors,
            "size": node_sizes,
            "opacity": node_alpha,
        },
        text=list(G.nodes()),
        hoverinfo="text+x+y+z",
    )

    # Loop on the list of edges to get the x,y,z, coordinates of the connected nodes
    # Those two points are the extrema of the line to be plotted
    x_edges = []
    y_edges = []
    z_edges = []

    for node_a, node_b in G.edges(data=False):
        x_edges.extend([pos[node_a][0], pos[node_b][0], None])
        y_edges.extend([pos[node_a][1], pos[node_b][1], None])
        z_edges.extend([pos[node_a][2], pos[node_b][2], None])

    axis = dict(
        showbackground=False,
        showline=False,
        zeroline=False,
        showgrid=False,
        showticklabels=False,
        title="",
    )

    repeated_edge_colours = []
    for (
        edge_col
    ) in (
        edge_colors
    ):  # Repeat as each line segment is ({x,y,z}_start, {x,y,z}_end, None)
        repeated_edge_colours.extend((edge_col, edge_col, edge_col))

    edge_colors = repeated_edge_colours

    edge_text = [
        " / ".join(list(edge_type))
        for edge_type in nx.get_edge_attributes(G, "kind").values()
    ]
    edge_text = np.repeat(
        edge_text, 3
    )  # Repeat as each line segment is ({x,y,z}_start, {x,y,z}_end, None)

    edges = go.Scatter3d(
        x=x_edges,
        y=y_edges,
        z=z_edges,
        mode="lines",
        line={"color": edge_colors, "width": 10},
        text=edge_text,
        hoverinfo="text",
    )

    return go.Figure(
        data=[nodes, edges],
        layout=go.Layout(
            title=plot_title,
            width=figsize[0],
            height=figsize[1],
            showlegend=False,
            scene=dict(
                xaxis=dict(axis),
                yaxis=dict(axis),
                zaxis=dict(axis),
            ),
            margin=dict(t=100),
        ),
    )


'''
Modified from graphein.protein.visualisation by Cam M
'''
def motif_plot_protein_structure_graph(
    G: nx.Graph,
    angle: int = 30,
    plot_title: Optional[str] = None,
    figsize: Tuple[int, int] = (7, 7),
    node_alpha: float = 0.7,
    node_size_min: float = 20.0,
    node_size_multiplier: float = 20.0,

    # modified by me:
    node_size_feature: str = "degree",
    # ---------------

    label_node_ids: bool = True,
    node_colour_map=plt.cm.plasma,
    edge_color_map=plt.cm.plasma,
    colour_nodes_by: str = "degree",
    colour_edges_by: str = "kind",
    edge_alpha: float = 0.5,
    plot_style: str = "ggplot",
    out_path: Optional[str] = None,
    out_format: str = ".png",
) -> Axes3D:
    """
    Plots protein structure graph in ``Axes3D``.

    :param G:  nx.Graph Protein Structure graph to plot.
    :type G: nx.Graph
    :param angle:  View angle. Defaults to ``30``.
    :type angle: int
    :param plot_title: Title of plot. Defaults to ``None``.
    :type plot_title: str, optional
    :param figsize: Size of figure, defaults to ``(10, 7)``.
    :type figsize: Tuple[int, int]
    :param node_alpha: Controls node transparency, defaults to ``0.7``.
    :type node_alpha: float
    :param node_size_min: Specifies node minimum size, defaults to ``20``.
    :type node_size_min: float
    :param node_size_multiplier: Scales node size by a constant. Node sizes reflect degree. Defaults to ``20``.
    :type node_size_multiplier: float
    :param node_size_feature: What feature to use to scale the node size. Defaults to ``degree``.
    :type node_size_feature: str
    :param label_node_ids: bool indicating whether or not to plot ``node_id`` labels. Defaults to ``True``.
    :type label_node_ids: bool
    :param node_colour_map: colour map to use for nodes. Defaults to ``plt.cm.plasma``.
    :type node_colour_map: plt.cm
    :param edge_color_map: colour map to use for edges. Defaults to ``plt.cm.plasma``.
    :type edge_color_map: plt.cm
    :param colour_nodes_by: Specifies how to colour nodes. ``"degree"``, ``"seq_position"`` or a node feature.
    :type colour_nodes_by: str
    :param colour_edges_by: Specifies how to colour edges. Currently only ``"kind"`` is supported.
    :type colour_edges_by: str
    :param edge_alpha: Controls edge transparency. Defaults to ``0.5``.
    :type edge_alpha: float
    :param plot_style: matplotlib style sheet to use. Defaults to ``"ggplot"``.
    :type plot_style: str
    :param out_path: If not none, writes plot to this location. Defaults to ``None`` (does not save).
    :type out_path: str, optional
    :param out_format: Fileformat to use for plot
    :type out_format: str
    :return: matplotlib Axes3D object.
    :rtype: Axes3D
    """

    # Get Node Attributes
    pos = nx.get_node_attributes(G, "coords")

    # Get node colours
    node_colors = colour_nodes(
        G, colour_map=node_colour_map, colour_by=colour_nodes_by
    )
    edge_colors = colour_edges(
        G, colour_map=edge_color_map, colour_by=colour_edges_by
    )

    # 3D network plot
    with plt.style.context(plot_style):

        fig = plt.figure(figsize=figsize)
        ax = Axes3D(fig, auto_add_to_figure=True)

        
        # TODO: incoroprate something like this dict:
        dict(degree=G.degree[key],
            asa=G.nodes[key],
        )

        # Get node scaling function
        def node_scale_size(G, feature):
            if feature == 'degree':
                return lambda k : node_size_min + node_size_multiplier * G.degree[k]
            elif feature in ['rsa', 'asa']:
                return lambda k : node_size_min + node_size_multiplier * G.nodes(data=True)[k]
            
        
        # Loop on the pos dictionary to extract the x,y,z coordinates of each node
        for i, (key, value) in enumerate(pos.items()):
            xi = value[0]
            yi = value[1]
            zi = value[2]

            # Scatter plot
            ax.scatter(
                xi,
                yi,
                zi,
                color=node_colors[i],
                s=node_scale_size(key),
                edgecolors="k",
                alpha=node_alpha,
            )
            if label_node_ids:
                label = list(G.nodes())[i]
                ax.text(xi, yi, zi, label)

        # Loop on the list of edges to get the x,y,z, coordinates of the connected nodes
        # Those two points are the extrema of the line to be plotted
        for i, j in enumerate(G.edges()):
            x = np.array((pos[j[0]][0], pos[j[1]][0]))
            y = np.array((pos[j[0]][1], pos[j[1]][1]))
            z = np.array((pos[j[0]][2], pos[j[1]][2]))

            # Plot the connecting lines
            ax.plot(x, y, z, c=edge_colors[i], alpha=edge_alpha)

    # Set title
    ax.set_title(plot_title)
    # Set the initial view
    ax.view_init(30, angle)
    # Hide the axes
    ax.set_axis_off()
    if out_path is not None:
        plt.savefig(out_path + str(angle).zfill(3) + out_format)
        plt.close("all")

    return ax

