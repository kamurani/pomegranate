#! /usr/bin/env python3
"""Load graphs from list of PDB files"""


"""
usage: python load.py ../datasets/yeast_full.txt structures/yeast-AF2 saved_graphs
"""

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



'''
    "-p",
    "--phosphosite",
    help="The file containing the phosphosite information",
    
)
@c.option(
    "--structures",
    help="The directory containing .PDB files.  Downloaded files will be stored here.",
    type=c.Path(
        exists=True, file_okay=False, dir_okay=True, path_type=pathlib.Path
    ),
) 

'''

@c.command()
@c.argument(
    'graphs', nargs=1,
)
@c.option(
    "-e",
    "--epochs",
    help="Number of epochs to train for",
    type=int,
    default=200,
)
#@c.argument('structures', nargs=1)
#@c.argument('graphs', nargs=1)
def main(
    graphs,
    epochs,
):

    graph_path = Path(graphs)

    

    # TODO: autodetect which file to use as infile, if directory is given.
    """
    if not graph_path.is_dir():
        raise ValueError(f"No such directory {graph_path}")
    
     
    filename = "graph_objects"
    in_path = os.path.join(graph_path, filename)

    in_path = Path(in_path)
    """

    if not graph_path.is_file():
        raise ValueError(f"Specified path {graphs} is not a file.")

    in_path = graph_path

    save_path = os.path.join(graph_path, "embeddings")

    print(f"Input file is {in_path}.")

    # Unpickle
    infile = open(in_path,'rb')
    loaded_graphs = pickle.load(infile)
    infile.close()

    loaded = list(loaded_graphs.values())

    graphs = []
    
    """
    for node_id, node_data in g.nodes(data=True):
        print (node_id)
        print (node_data)
    return
    """
    for i in range(len(loaded)):
        g = loaded[i]['graph'].copy()
        #print(g.nodes())
        psite = loaded[i]['psite']
        res = loaded[i]['res']
        #print(res)
        #print("psite coords:")
        #print(psite['coords'])
        #print(res in g.nodes())

        psite_loc = np.array(psite['coords'])



        #print(psite_coords)
        
        # TODO: use EDGE FEATURES IN GRAPH LOADING TOO
        
        # For each node, create `feature` vector. 
        for node_id, node_data in g.nodes(data=True):
            # Create list of numerical features for each node.
            #node_data = g.nodes(node_id)
            
            m = node_data["meiler"] 
            node_loc = np.array(node_data["coords"])
            dist_to_psite = np.linalg.norm(psite_loc - node_loc)
            #print(f"dist to psite: {dist_to_psite}")

            rsa = node_data["rsa"]
            b_fac = node_data["b_factor"]
            #feature = [*m, *c, rsa]

            feature = [dist_to_psite, rsa, b_fac, *m]  # 10 long   

            node_data["feature"] = feature
            
        # Create sg instance from graph, using `feature` vector. 
        g_attr = StellarGraph.from_networkx(g, node_features="feature")
        graphs.append(g_attr)

    
    #print(graphs[5].info())

    # Generator
    generator = sg.mapper.PaddedGraphGenerator(graphs)

    gc_model = sg.layer.GCNSupervisedGraphClassification(
        [8, 4], ["relu", "relu"], generator, pool_all_layers=True
    )
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

    train_gen = generator.flow(graph_idx, batch_size=16, targets=targets)

    pair_model.compile(keras.optimizers.Adam(1e-2), loss="mse")

    NUM_EPOCHS = epochs # 500
    history = pair_model.fit(train_gen, epochs=NUM_EPOCHS, verbose=0)
    #sg.utils.plot_history(history)

    embeddings = embedding_model.predict(generator.flow(graphs))

    print(f"Embeddings have shape {embeddings.shape}")

    # Save embeddings 
    print("Saving embeddings...")
    outfile =  open(save_path, 'wb')
    pickle.dump(embeddings, outfile)
    outfile.close()
    print(f"Saved embeddings at {outfile}.")




    
    # Plot 
    from sklearn.manifold import TSNE

    tsne = TSNE(2, init='random')
    two_d = tsne.fit_transform(embeddings)

    from matplotlib import pyplot as plt

    plt.scatter(two_d[:, 0], two_d[:, 1], alpha=0.6)

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