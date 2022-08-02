#! /usr/bin/env python3
"""Cluster embeddings"""

import pickle
from protein.phosphosite import get_surface_motif


import pathlib
from graphein.protein.graphs import construct_graph
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.edges.distance import *	

from graphein.protein.config import DSSPConfig
from graphein.protein.features.nodes import rsa

from pathlib import Path

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

"""
@c.option(
    "-e",
    "--epochs",
    help="Number of epochs to train for",
    type=int,
    default=200,
)
#@c.argument('structures', nargs=1)
#@c.argument('graphs', nargs=1)
"""


@c.command()
@c.argument(
    'embeddings', nargs=1,
)
def main(
    embeddings,
):

    em_path = Path(embeddings)
    if not em_path.is_dir():
        raise ValueError(f"No such directory {em_path}")
    
     
    filename = "embeddings"
    in_path = os.path.join(em_path, filename)

    save_path = os.path.join(em_path, "embedding_clustered")

    print(f"Input file is {in_path}.")

    # Unpickle
    infile = open(in_path,'rb')
    embeddings = pickle.load(infile)
    infile.close()

    print("EMBEDDINGS:")
    print(embeddings.shape)

    print(embeddings[0])

    return

    # Save embeddings 
    '''
    print("Saving embeddings...")
    outfile =  open(save_path, 'wb')
    pickle.dump(embeddings, outfile)
    outfile.close()
    print(f"Saved embeddings at {outfile}.")

    '''

    

    
    # Plot 

    from sklearn.manifold import TSNE

    tsne = TSNE(2, init='random')
    two_d = tsne.fit_transform(embeddings)

    from matplotlib import pyplot as plt

    plt.scatter(two_d[:, 0], two_d[:, 1], alpha=0.6)

    






if __name__ == "__main__":
    main()