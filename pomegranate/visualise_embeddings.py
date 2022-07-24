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
@c.argument(
    'outdir', nargs=1,

)
@c.option(
    "-N",
    "--num-plots",
    help="How many plots to output.",
    type=c.INT,
    default=1, 
    show_default=True,
)
@c.option(
    # TODO: support multiple ? i.e. PCA first?
    "--dim-method",
    "--method",
    type=c.Choice(['UMAP', 'tSNE', 'PCA'], case_sensitive=False),
    help="Method to use for dimensionality reduction.",
    default="tSNE",
    show_default=True,
)
@c.option(
    "--labels/--no-labels",
    "-l",
    default=True,
    help="Show labels on plot."
)
@c.option(
    "-i",
    "--interactive",
    is_flag=True,
    help="Open interactive plot viewer.",
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
    outdir,

    # OPTIONAL ARGUMENTS
    num_plots,

    # OPTIONS
    dim_method,
    labels,
    interactive,
    verbose,
):
    verbose = True # TODO: remove

    label_graph = labels # rename variable

    in_path = Path(embeddings)
    out_dir = Path(outdir)

    

    # TODO: autodetect which file to use as infile, if directory is given.
    
    if not out_dir.is_dir():
        raise ValueError(f"No such directory {out_dir}")
    
    """
    filename = "graph_objects"
    in_path = os.path.join(graph_path, filename)

    in_path = Path(in_path)
    """

    if not in_path.is_file():
        raise ValueError(f"Specified path {embeddings} is not a file.")


    basename = os.path.basename(in_path)
    parent = os.path.dirname(in_path)
    save_path = os.path.join(out_dir, f"{basename}_PLOT_X.png")

    print(f"Input file is {in_path}.  Saving plots to {save_path}")

    

    # Unpickle
    if verbose: print(f"Loading...", end=" ")
    infile = open(in_path,'rb')
    data = pickle.load(infile)
    infile.close()

    embeddings  = data['embeddings']
    labels      = data['labels']


    if verbose: print(f"Loaded embeddings of shape {embeddings.shape}, labels of length {len(labels)}")

    """
    for node_id, node_data in g.nodes(data=True):
        print (node_id)
        print (node_data)
    return
    """

    
    kinases = np.array(labels)
    kinases = np.unique(kinases)    # unique kinases

    # get dict mapping kinase to number
    mapping = {}
    for i, k in enumerate(kinases):
        mapping[k] = i

    group = []
    for i, l in enumerate(labels):
        group.append(mapping[l])
    
    # Plot 
    if dim_method.lower() == "tsne":
        from sklearn.manifold import TSNE
        tsne = TSNE(
            2, 
            init='pca', # 'random'
            learning_rate='auto'
        )
        two_d = tsne.fit_transform(embeddings)
    else:
        raise NotImplementedError(f"Dimensionality reduction using '{dim_method}' not implemented.")

    from matplotlib import pyplot as plt



    for i in range(num_plots):

        fig, ax = plt.subplots()

        # TSNE random
        tsne = TSNE(
            2, 
            init='random',
            learning_rate='auto'
        )
        two_d = tsne.fit_transform(embeddings)

        plt.scatter(
            two_d[:, 0], 
            two_d[:, 1],
            c=group,
            label=group,
            alpha=0.6
        )

        if label_graph:
            for i, l in enumerate(labels):
                plt.text(
                    x=two_d[i,0],
                    y=two_d[i,1],
                    s=l,
                )

        #for i, txt in enumerate(graph_labels):
        #    ax.annotate(txt, two_d[i, 0], two_d[i, 1])

        save_path = os.path.join(out_dir, f"{basename}_PLOT_{i}.png")
        if verbose: print(f"[{i}] Saving...", end=" ")
        try:
            plt.savefig(save_path)
            if verbose: print(f"DONE at {save_path}")
        except:
            if verbose: print(f"FAILED.")

    if interactive:
        plt.show()
    
    return

    
    

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