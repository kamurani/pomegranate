"""Custom plot functions.  Modified from graphein.protein.visualisation"""


from __future__ import annotations

import logging
from itertools import count
from typing import Dict, List, Optional, Tuple, Union

import plotly
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.distance import pdist, squareform
import re

from graphein.protein.subgraphs import extract_k_hop_subgraph
from graphein.utils.utils import import_message
import graphein.protein.edges.distance as g_dist

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
    if g is None:
        raise ValueError("Must provide a graph as input.")

    # Get distance matrix
    dist_mat = ordered_distmat(g, aa_order)      
    
    # Get labels
    x_range = list(dist_mat.index)
    y_range = list(dist_mat.columns)
    
    # add phospho label to selected site
    for i in range(len(x_range)):
        s = x_range[i]
        if psite == s:
            break
    
    x_range[i] = f"{x_range[i]} [P]"
    y_range[i] = f"[P] {y_range[i]}"
    
    # Default title
    if not title:
        title = g.graph["name"] + " - Distance Matrix"


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



# TODO: modify which attributes we use in node's label
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
    node_size_feature: str = "degree",
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
    :param node_size_feature: What feature to use to scale the node size by. Defaults to ``degree``.
    :type node_size_feature: str
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

    # Get node size 
    def node_scale_by(G, feature):
        if feature == 'degree':
            return lambda k : node_size_min + node_size_multiplier * G.degree[k]
        elif feature == 'rsa':
            return lambda k : node_size_min + node_size_multiplier * G.nodes(data=True)[k]['rsa']
        # Meiler embedding dimension
        p = re.compile("meiler-([1-7])")
        dim = p.search(feature).group(1)
        if dim:
            return lambda k : node_size_min + node_size_multiplier * max(0, G.nodes(data=True)[k]['meiler'][f'dim_{dim}']) # Meiler values may be negative
        else:
            raise ValueError(f"Cannot size nodes by feature '{feature}'")   

    get_node_size = node_scale_by(G, node_size_feature) 

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
        node_sizes.append(get_node_size(key))

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
        def node_scale_by(G, feature):
            if feature == 'degree':
                return lambda k : node_size_min + node_size_multiplier * G.degree[k]
            elif feature in ['rsa', 'asa']:
                return lambda k : node_size_min + node_size_multiplier * G.nodes(data=True)[k]
        
        node_scale_size = node_scale_by(G, node_size_feature)
        
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



"""
Modified from graphein.protein.visualisation by Cam M
"""

def motif_asteroid_plot(
    g: nx.Graph,
    node_id: str,
    k: int = 2,
    colour_nodes_by: str = "shell",  # residue_name, hydrophobicity
    size_nodes_by: str = "degree",   # RSA 
    colour_edges_by: str = "kind",
    edge_colour_map: plt.cm.Colormap = plt.cm.plasma,
    edge_alpha: float = 1.0,
    show_labels: bool = True,
    title: Optional[str] = None,
    width: int = 600,
    height: int = 500,
    use_plotly: bool = True,
    show_edges: bool = False,
    show_legend: bool = True,
    node_size_min: float = 20,
    node_size_multiplier: float = 10,
) -> Union[plotly.graph_objects.Figure, matplotlib.figure.Figure]:
    """Plots a k-hop subgraph around a node as concentric shells.

    Radius of each point is proportional to a node attribute, with degree as default. (modified by node_size_multiplier).

    :param g: NetworkX graph to plot.
    :type g: nx.Graph
    :param node_id: Node to centre the plot around.
    :type node_id: str
    :param k: Number of hops to plot. Defaults to ``2``.
    :type k: int
    :param colour_nodes_by: Colour the nodes by this attribute. 
    :type colour_nodes_by: str


    
    :param size_nodes_by: Size the nodes by this attribute.  
    :type size_nodes_by: str

    :param colour_edges_by: Colour the edges by this attribute. Currently only ``"kind"`` is supported.
    :type colour_edges_by: str
    :param edge_colour_map: Colour map for edges. Defaults to ``plt.cm.plasma``.
    :type edge_colour_map: plt.cm.Colormap
    :param edge_alpha: Sets a given alpha value between 0.0 and 1.0 for all the edge colours.
    :type edge_alpha: float
    :param title: Title of the plot. Defaults to ``None``.
    :type title: str
    :param width: Width of the plot. Defaults to ``600``.
    :type width: int
    :param height: Height of the plot. Defaults to ``500``.
    :type height: int
    :param use_plotly: Use plotly to render the graph. Defaults to ``True``.
    :type use_plotly: bool
    :param show_edges: Whether to show edges in the plot. Defaults to ``False``.
    :type show_edges: bool
    :param show_legend: Whether to show the legend of the edges. Fefaults to `True``.
    :type show_legend: bool
    :param node_size_min: Specifies node minimum size. Defaults to ``20.0``.
    :type node_size_min: float
    :param node_size_multiplier: Multiplier for the size of the nodes. Defaults to ``10``.
    :type node_size_multiplier: float.
    :returns: Plotly figure or matplotlib figure.
    :rtpye: Union[plotly.graph_objects.Figure, matplotlib.figure.Figure]
    """
    assert node_id in g.nodes(), f"Node {node_id} not in graph"

    nodes: Dict[int, List[str]] = {0: [node_id]}
    node_list: List[str] = [node_id]

    # Iterate over the number of hops and extract nodes in each shell
    for i in range(1, k):
        subgraph = extract_k_hop_subgraph(g, node_id, k=i)
        candidate_nodes = subgraph.nodes()
        # Check we've not already found nodes in the previous shells
        nodes[i] = [n for n in candidate_nodes if n not in node_list]
        node_list += candidate_nodes
    shells = [nodes[i] for i in range(k)]


    #log.debug(f"Plotting shells: {shells}")

    if use_plotly:
        # Get shell layout and set as node attributes.
        pos = nx.shell_layout(subgraph, shells)
        nx.set_node_attributes(subgraph, pos, "pos")

        if show_edges:
            edge_colors = colour_edges(
                subgraph,
                colour_map=edge_colour_map,
                colour_by=colour_edges_by,
                set_alpha=edge_alpha,
                return_as_rgba=True,
            )
            show_legend_bools = [
                (True if x not in edge_colors[:i] else False)
                for i, x in enumerate(edge_colors)
            ]
            edge_trace = []
            for i, (u, v) in enumerate(subgraph.edges()):
                x0, y0 = subgraph.nodes[u]["pos"]
                x1, y1 = subgraph.nodes[v]["pos"]
                bond_kind = " / ".join(list(subgraph[u][v]["kind"]))
                tr = go.Scatter(
                    x=(x0, x1),
                    y=(y0, y1),
                    mode="lines",
                    line=dict(width=1, color=edge_colors[i]),
                    hoverinfo="text",
                    text=[bond_kind],
                    name=bond_kind,
                    legendgroup=bond_kind,
                    showlegend=show_legend_bools[i],
                )
                edge_trace.append(tr)

        node_x: List[str] = []
        node_y: List[str] = []
        for node in subgraph.nodes():
            x, y = subgraph.nodes[node]["pos"]
            node_x.append(x)
            node_y.append(y)


        
        
        def node_size_function(g: nx.Graph, feature: str):
            if feature == 'degree':
                return lambda k : g.degree(k)
            elif feature == 'rsa':
                return lambda k : g.nodes(data=True)[k]['rsa']
            else:
                raise NotImplementedError(
                    f"Size by {size_nodes_by} not implemented."
                )

        node_size = node_size_function(subgraph, size_nodes_by)

        node_sizes = [
            node_size_min + node_size(n) * node_size_multiplier for n in subgraph.nodes()
        ]

        '''
        if size_nodes_by == "degree":
            node_sizes = [
                node_size_min + subgraph.degree(n) * node_size_multiplier for n in subgraph.nodes()
            ]
        elif size_nodes_by == "rsa":
            # TODO: check that graph actually has rsa attribute
            print("Sizing nodes by rsa...")
            node_sizes = [
                subgraph.nodes(data=True)[n]['rsa'] * node_size_multiplier for n in subgraph.nodes()
            ]
            
        else:
            raise NotImplementedError(
                f"Size by {size_nodes_by} not implemented."
            )
        '''

        if colour_nodes_by == "shell":
            node_colours = []
            for n in subgraph.nodes():
                for k, v in nodes.items():
                    if n in v:
                        node_colours.append(k)
                        print(f"k: {k}")

        # IMPLEMENTED BY ME:
        elif colour_nodes_by == "hydrophobicity":
            node_colours = []
            for n in subgraph.nodes():
                for k, v in nodes.items():
                    if n in v:
                        print(f"n: {n}")
                        node_colours.append(aa2hydrophobicity(n.split(':')[1]))
        
        # TODO
        elif colour_nodes_by == "residue_name":
            node_colours = []
            for n in subgraph.nodes():
                pass
        else:
            raise NotImplementedError(
                f"Colour by {colour_nodes_by} not implemented."
            )
            # TODO colour by AA type
        node_trace = go.Scatter(
            x=node_x,
            y=node_y,
            text=list(subgraph.nodes()),
            mode="markers+text" if show_labels else "markers",
            hoverinfo="text",
            textposition="bottom center",
            showlegend=False,
            marker=dict(
                colorscale="viridis",
                reversescale=True,
                color=node_colours,
                size=node_sizes, 
                colorbar=dict(
                    thickness=15,
                    title=str.capitalize(colour_nodes_by),
                    tickvals=list(range(k)),
                    xanchor="left",
                    titleside="right",
                ),
                line_width=2,
            ),
        )

        data = edge_trace + [node_trace] if show_edges else [node_trace]
        fig = go.Figure(
            data=data,
            layout=go.Layout(
                title=title if title else f'Asteroid Plot - {g.graph["name"]}',
                width=width,
                height=height,
                titlefont_size=16,
                legend=dict(yanchor="top", y=1, xanchor="left", x=1.10),
                showlegend=True if show_legend else False,
                hovermode="closest",
                margin=dict(b=20, l=5, r=5, t=40),
                xaxis=dict(
                    showgrid=False, zeroline=False, showticklabels=False
                ),
                yaxis=dict(
                    showgrid=False, zeroline=False, showticklabels=False
                ),
            ),
        )
        return fig
    else:
        nx.draw_shell(subgraph, nlist=shells, with_labels=show_labels)


'''
Distance matrix in other orders
Modified from graphein.protein.edges.distance by Naomi Warren
'''
def ordered_distmat(g: nx.graph, order: str) -> pd.DataFrame:
    """
    Compute pairwise Euclidean distances between every atom, ordering the matrix
    by 'order'.
    :raises: ValueError if ``g.graph['pdb_df']`` does not contain the required columns.
    :return: pd.Dataframe of Euclidean distance matrix.
    :rtype: pd.DataFrame
    """

    pdb_df = g.graph['pdb_df']

    # Check df has the correct columns (not sure why it wouldn't? but ...)
    if (
        not pd.Series(["x_coord", "y_coord", "z_coord"])
        .isin(pdb_df.columns)
        .all()
    ):
        raise ValueError(
            "Dataframe must contain columns ['x_coord', 'y_coord', 'z_coord']"
        )

    # Get the distances as a square matrix
    eucl_dists = pdist(
        pdb_df[["x_coord", "y_coord", "z_coord"]], metric="euclidean"
    )
    eucl_dists = pd.DataFrame(squareform(eucl_dists))

    # Re-order the rows in the appropriate order
    eucl_dists.index = pdb_df.node_id
    eucl_dists.columns = pdb_df.node_id

    cur_order = list(g.nodes)
    # TODO: seems to be issue where g.nodes isn't consistent with g.graph['pdb_df'].
    # TODO: need to work out why/where one is being modified and not the other.
    print (f'Length of pdb_df: {len(pdb_df)}')
    print (f'Length of g.nodes: {len(cur_order)}')

    if order == 'seq':
        # No changes required
        return eucl_dists
    elif order == 'hydro':
        # Get the order from low to high hydrophobicity
        ordering = "IVLFCMAWGTSYPHNDQEKR"
        hydro_sorted = sorted(cur_order, key=lambda res: [ordering.index(a) for a in aa3to1(res.split(':')[1])])
        # Re-index (i.e. switch rows and columns) in this order
        eucl_dists = eucl_dists.reindex(hydro_sorted)
        eucl_dists = eucl_dists.reindex(columns=hydro_sorted)

    return eucl_dists

"""
Get hydrophobicity map
"""
def aa2hydrophobicity(
    aa: str,
    mapping : str = 'a',
    ):
    if mapping == 'a':
        hmap = { 
            "ILE" : 4.5,
            "VAL" : 4.2,
            "LEU" : 3.8,
            "PHE" : 2.8,
            "CYS" : 2.5,
            "MET" : 1.9,
            "ALA" : 1.8,
            "GLY" : -0.4,
            "THR" : -0.7,
            "SER" : -0.8,
            "TRP" : -0.9,
            "TYR" : -1.3,
            "PRO" : -1.6,
            "HIS" : -3.2,
            "GLU" : -3.5,
            "GLN" : -3.5,
            "ASP" : -3.5,
            "ASN" : -3.5,
            "LYS" : -3.9,
            "ARG" : -4.5,
        }
    return hmap[aa]

