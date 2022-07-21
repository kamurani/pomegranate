#! /usr/bin/env python3
"""Load graphs from list of PDB files"""


"""
usage: python load.py ../datasets/yeast_full.txt structures/yeast-AF2 saved_graphs
"""

import collections
import inspect
import pickle
import re
import time
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

# TODO: support for describing features like 'm1-4' instead of 'm1','m2'...
GRAPH_NODE_FEATURES = [
    'm1','m2','m3','m4','m5','m6','m7'
    'pdist', 'coords', 'bfac', 'rsa',
]
GRAPH_EDGE_FEATURES = []





"""
Flatten an irregular list of iterables / ints. 
"""
def flatten(xs):
    for x in xs:
        if isinstance(x, Iterable) and not isinstance(x, (str, bytes)):
            yield from flatten(x)
        else:
            yield x

def nx_to_sg(
    graphs: Dict[int, Dict[str, Union[str, nx.Graph]]] = None,
    include_features: List[str] = None,
    verbose: bool = False,
) -> Dict[int, Dict[str, Union[str, sg.StellarGraph]]]:
    in_graphs: List[Dict[str, Union[str, nx.Graph]]] = list(graphs.values())

    out_graphs = {}

    for idx, graph in enumerate(in_graphs):

        if graph is None: continue
        if graph['graph'] is None: continue
        g       = graph['graph'].copy()
        kinase  = graph['kinase']
        psite   = graph['psite']
        res     = graph['res']

        
        

        # GET LOCATION OF PSITE
        psite_coords = np.array(psite['coords'])

        # TODO: 'zero' coords to use psite as (0,0,0) and all other nodes are 'relative' vectors from here
        # TODO: use EDGE FEATURES IN GRAPH LOADING TOO
        graph_labels = []

        """
        dispatcher = {'foobar': [foo, bar], 'bazcat': [baz, cat]}
        
        def fire_all(func_list):
            for f in func_list:
                f()
        
        fire_all(dispatcher['foobar])
        """

        # WHAT IF?
        '''
        dispatcher = {
            'rsa' : get_feature_func('rsa'), 
            'm' : get_feature_func('meiler'),
            'coords', get_coords,
        }
        '''
        dispatcher: Dict[str, Callable] = {}

        def get_rsa(node) -> float:
            return node["rsa"] 
        def get_bfac(node) -> float:
            return node["b_factor"]
        def get_dist_psite(node) -> float:
            return 0 # TODO
        def get_coords(node) -> np.ndarray:
            return np.array(node["coords"])   # np array
        def get_meiler_func(m: str) -> Callable:
            p = re.compile("m([1-7])")
            dim = p.search(m).group(1)
            if dim:
                idx = int(dim) - 1
                return lambda node_data : node_data["meiler"][idx]
            else: 
                raise ValueError(f"No such meiler embedding dimension '{m}'")
        
        '''
        Returns a function that calculates distance to a phosphosite for any node that is input.
        '''
        def get_psite_func(psite_coords: np.ndarray) -> Callable:
            return lambda node_data: np.linalg.norm(get_coords(node_data) - psite_coords)
        
        dispatcher = {
            'rsa' : get_rsa, 
            'coords': get_coords,
            'bfac': get_bfac, 
            'pdist': get_psite_func(psite_coords),
        }
        disp = dispatcher

        for m in ['m1', 'm2', 'm3', 'm4', 'm5', 'm6', 'm7']: # each meiler embedding dimension
            disp[m] = get_meiler_func(m)

        
        if include_features is None:   
            # Use default features     
            include_features = ['rsa', 'pdist', 'bfac', 'm1', 'm2', 'm3']
            if verbose:
                print(f"No feature list given.  Using default: {include_features}")

        # Check features given
        for f in include_features: assert f in disp.keys(), f"Unknown node feature {f}"

        for node_id, node_data in g.nodes(data=True):
            
            n = node_data
            if True:
                
                n['feature'] = list(flatten([disp[f](n) for f in include_features]))

                print(f"Feature: {n['feature']}")
                
                # WHY DOESN'T THIS WORK? TODO
                #node_data['feature'] = [*disp[f](n) if is_iterable(disp[f](n)) else disp[f](n) for f in features]

            else:

                # Get list of numerical features
                m = node_data["meiler"] # node's meiler embedding
                node_loc = np.array(node_data["coords"])    # node's euclidean coords
                dist_to_psite = np.linalg.norm(node_loc - psite_coords) # node's distance to psite
                
                rsa = node_data["rsa"]  # node's relative solvent accessibility
                b_fac = node_data["b_factor"]   # node's temperature 
                
                #feature = [*m, *c, rsa]
                #feature = [dist_to_psite, rsa, b_fac, *m]  # 10 long   
                
                feature = [dist_to_psite, rsa, b_fac, m[2], m[3], m[4], m[5]]  # 10 long 
                #node_data["feature"] = feature  # create new 'feature' vector to be accessed by StellarGraph instance
            
        # Create sg instance from graph, using `feature` vector. 
        g_attr = StellarGraph.from_networkx(g, node_features="feature")

        # Create 'graph' dict
        out_graphs[idx] = dict(graph=g_attr, res=res, psite=psite, kinase=kinase)
        
        graph_labels.append(kinase) # unused for now

    return out_graphs

"""
Get just the 'graph' components of a graph dict.
"""
def get_graphs_and_labels(
    graphs: Dict[int, Dict[str, Union[str, sg.StellarGraph, nx.Graph]]]
) -> Tuple[List[Union[sg.StellarGraph, nx.Graph]], List[str]]:

    labels          = []
    graph_objects   = []
    graph_dicts     = []
    
    for graph in graphs.values():

        if isinstance(graph['graph'], (sg.StellarGraph, nx.Graph)):
            graph_objects.append(graph['graph'])
            graph_dicts.append(graph)
            labels.append(graph['kinase'])
        

    return graph_objects, graph_dicts, labels

"""
Get list of graph labels from a graph dict object
"""
def get_graph_labels(
    graphs: Dict[int, Dict[str, Union[str, sg.StellarGraph, nx.Graph]]]
) -> List[str]:

    labels = []

    for i in range(len(graphs)):    
        if graphs[i] is None:
            labels.append(None) # preserve indexes, so we append even if None
        else:
            labels.append(graphs[i]['kinase'])
      

"""
Test that ``graph_labels`` array correctly describes a ``graphs`` dict
"""
def test_labels(graphs: Dict, labels: List):
    for i in graphs.keys():
        assert graphs[i]['kinase'] == labels[i], "Looks like we screwed up graph labels... yikes"

"""
Test ``graph_labels`` array correctly describes a list of ``graph_dict``
"""
def test_labels_list(graphs: List, labels: List):
    for i, graph in enumerate(graphs):
        assert graph['kinase'] == labels[i]


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
    help="Number of epochs to train for.",
    type=c.INT,
    default=200,
    show_default=True,
)
@c.option(
    "-b",
    "--batch-size",
    type=c.INT,
    default=16, show_default=True,

    help="Batch size to be used by data generator."
)
@c.option(
    "-v", 
    "--verbose",
    is_flag=True,
    
)
#@c.argument('structures', nargs=1)
#@c.argument('graphs', nargs=1)
def main(
    graphs,
    epochs,
    batch_size,
    verbose,
):
    verbose = True # TODO: remove


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

    parent = os.path.dirname(in_path)
    save_path = os.path.join(parent, "embeddings")

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

    graphs, graph_dicts, graph_labels = get_graphs_and_labels(graphs_dict)


    test_labels_list(graph_dicts, graph_labels)
    
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

    num_samples = 600
    graph_idx = np.random.RandomState(0).randint(len(graphs), size=(num_samples, 2))

    targets = [graph_distance(graphs[left], graphs[right]) for left, right in graph_idx]

    train_gen = generator.flow(graph_idx, batch_size=batch_size, targets=targets)

    pair_model.compile(keras.optimizers.Adam(1e-2), loss="mse")
    
    if verbose:
        print("-----------------")
        print(f"Starting training on {num_samples} samples for {epochs} epochs...")
    
    start = time.time()
    history = pair_model.fit(train_gen, epochs=epochs, verbose=0) # verbose?
    end = time.time()
    #sg.utils.plot_history(history)

    if verbose:
        print("-----------------")
        print(f"Completed training after {end - start} seconds.\n")

    if verbose: print(f"Generating embeddings...", end=" ")
    embeddings = embedding_model.predict(generator.flow(graphs))
    if verbose: print(f"DONE.")
    print(f"Embeddings have shape {embeddings.shape}")

    # Save embeddings 
    data = dict(
        embeddings=embeddings,
        labels=graph_labels
    )
    if verbose: print("Saving embeddings...", end=" ")
    outfile =  open(save_path, 'wb')
    pickle.dump(data, outfile)
    outfile.close()
    if verbose: print(f"DONE.")
    print(f"Saved embeddings at {save_path}.")

    return

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