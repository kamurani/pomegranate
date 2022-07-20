#! /usr/bin/env python3
"""Visualise saved embeddings"""


"""
usage: python visualise_embeddings.py savedgraphs/embeddings 
"""

import collections
import inspect
import pickle
import re
from protein.phosphosite import get_surface_motif


import pathlib
from graphein.protein.graphs import construct_graph
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.edges.distance import *	

from graphein.protein.config import DSSPConfig
from graphein.protein.features.nodes import rsa

from pathlib import Path
from typing import Callable, List, Dict, Union, Tuple
from collections.abc import Iterable

import click as c




import stellargraph as sg
from stellargraph import StellarGraph

import pandas as pd
import numpy as np
import networkx as nx

from networkx import Graph
import tensorflow as tf
from tensorflow import keras

from IPython.display import display, HTML


import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'


import os
import urllib as ul

from definitions import STRUCTURE_PATH
from traitlets import default
from validate import get_database

# TODO: alternative AA residue property feature vectors
# e.g. 1-hot encoding of Residue type


# TODO: pass in directory to save the plots. 
# TODO: multiple plots in one go i.e. 3x3 grid to give tSNE 9 chances.
@c.command()
@c.argument(
    'embeddings', nargs=1,

)
@c.option(
    # TODO: support multiple ? i.e. PCA first?
    "--dim-method",
    "--method"
    type=c.Choice(['UMAP', 'tSNE', 'PCA'], case_sensitive=False),
    help="Method to use for dimensionality reduction.",
    default="tSNE",
    show_default=True,
)
@c.option(
    "--labels",
    "-l",
    is_flag=True,
    help="Show labels on plot."
)
@c.option(
    "-v", 
    "--verbose",
    is_flag=True,
)
#@c.argument('structures', nargs=1)
#@c.argument('graphs', nargs=1)
def main(
    # POSITIONAL ARGUMENTS
    embeddings,

    # OPTIONAL ARGUMENTS

    # OPTIONS
    dim_method,
    labels,
    verbose,
):
    verbose = True # TODO: remove


    in_path = Path(embeddings)

    

    # TODO: autodetect which file to use as infile, if directory is given.
    """
    if not graph_path.is_dir():
        raise ValueError(f"No such directory {graph_path}")
    
     
    filename = "graph_objects"
    in_path = os.path.join(graph_path, filename)

    in_path = Path(in_path)
    """

    if not in_path.is_file():
        raise ValueError(f"Specified path {embeddings} is not a file.")

    parent = os.path.dirname(in_path)
    save_path = os.path.join(parent, "vis_plot.png")

    print(f"Input file is {in_path}.")

    # Unpickle
    infile = open(in_path,'rb')
    data = pickle.load(infile)
    graphs_dict = data['graphs_dict']
    graph_format = data['format']
    infile.close()

    loaded = list(graphs_dict.values())

    # Check which type of graph we are using. 
    g = loaded[0]['graph']
    if isinstance(g, sg.StellarGraph):
        print(f"StellarGraph format.")
    
    elif isinstance(g, nx.Graph):
        print(f"NetworkX format.")
        print(f"Converting...")

        raise NotImplementedError(f"Converting at this point not implemented yet.")

    
    """
    for node_id, node_data in g.nodes(data=True):
        print (node_id)
        print (node_data)
    return
    """

    # Get graphs and labels as lists
    graphs: List[sg.StellarGraph]
    graph_labels: List[str]

    graphs, graph_labels = get_graphs_and_labels(graphs_dict)
    
    # Generator
    generator = sg.mapper.PaddedGraphGenerator(graphs)

    # Model
    layers = [32, 16]
    act = "relu"
    activations = [act for i in range(len(layers))]
    gc_model = sg.layer.GCNSupervisedGraphClassification(
        #[32, 16, 8], ["relu", "relu", "relu"], generator, pool_all_layers=True
        layers, activations, generator, pool_all_layers=True
    )

    if verbose: 
        print(f"MODEL\n-----")
        print(f"LAYERS:\t\t{layers}")
        print(f"ACTIVATIONS:\t{activations}")

    inp1, out1 = gc_model.in_out_tensors()
    inp2, out2 = gc_model.in_out_tensors()

    vec_distance = tf.norm(out1 - out2, axis=1)

    pair_model = keras.Model(inp1 + inp2, vec_distance)
    embedding_model = keras.Model(inp1, out1)

    def graph_distance(graph1, graph2):
        spec1 = nx.laplacian_spectrum(graph1.to_networkx(feature_attr="feature"))
        spec2 = nx.laplacian_spectrum(graph2.to_networkx(feature_attr="feature"))
        k = min(len(spec1), len(spec2))
        dist = np.linalg.norm(spec1[:k] - spec2[:k])
        #print(f"Dist: {dist}")
        return dist 

    graph_idx = np.random.RandomState(0).randint(len(graphs), size=(200, 2))

    targets = [graph_distance(graphs[left], graphs[right]) for left, right in graph_idx]

    train_gen = generator.flow(graph_idx, batch_size=batch_size, targets=targets)

    pair_model.compile(keras.optimizers.Adam(1e-2), loss="mse")

    history = pair_model.fit(train_gen, epochs=epochs, verbose=0)
    #sg.utils.plot_history(history)

    if verbose: print(f"Generating embeddings...", end=" ")
    embeddings = embedding_model.predict(generator.flow(graphs))
    if verbose: print(f"DONE.")
    print(f"Embeddings have shape {embeddings.shape}")

    # Save embeddings 
    if verbose: print("Saving embeddings...", end=" ")
    outfile =  open(save_path, 'wb')
    pickle.dump(embeddings, outfile)
    outfile.close()
    if verbose: print(f"DONE.")
    print(f"Saved embeddings at {save_path}.")


    kinases = np.array(graph_labels)
    kinases = np.unique(kinases)    # unique kinases

    # get dict mapping kinase to number
    mapping = {}
    for i, k in enumerate(kinases):
        mapping[k] = i

    group = []
    for i, l in enumerate(graph_labels):
        group.append(mapping[l])
    
    # Plot 
    from sklearn.manifold import TSNE

    tsne = TSNE(2, learning_rate='auto')
    two_d = tsne.fit_transform(embeddings)

    from matplotlib import pyplot as plt

    fig, ax = plt.subplots()

    

    plt.scatter(
        two_d[:, 0], 
        two_d[:, 1],
        c=group,
        label=group,
        alpha=0.6
    )


    label_graph = True
    if label_graph:
        for i, l in enumerate(graph_labels):
            plt.text(
                x=two_d[i,0],
                y=two_d[i,1],
                s=l,
            )



    #for i, txt in enumerate(graph_labels):
    #    ax.annotate(txt, two_d[i, 0], two_d[i, 1])

    plt.show()
    plt.savefig('EMBEDDINGS_PLOT.png')


    
    

    '''

    # Convert to sg instances
    print("Converting...")
    for k in loaded_graphs.keys():
        
        g = loaded_graphs[k].copy()
        #print(f"[{k}] {g}")

        print(g)

        for node_id, node_data in g.nodes(data=True):
            # Create list of numerical features for each node.
            m = node_data["meiler"]
            c = node_data["coords"]
            rsa = node_data["rsa"]
            feature = [*m, *c, rsa]

            node_data["feature"] = feature
        
        g_attr = StellarGraph.from_networkx(g, node_features="feature")

        loaded_graphs[k] = StellarGraph.from_networkx(g_attr)
        #print(g.info())

    graphs = list(loaded_graphs.values())
    

    def compute_features(node_id):
    # in general this could compute something based on other features, but for this example,
    # we don't have any other features, so we'll just do something basic with the node_id
        hydromap_a = { 
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
        seq_pos = node_id.split(':')[-1]
        res = node_id.split(':')[1]
        hydro = hydromap_a[res] 
        return [seq_pos, hydro]

    '''

    


    
    
    

   
    



    




if __name__ == "__main__":
    main()