#! /usr/bin/env python3
"""Load graphs from list of PDB files"""


"""
usage: python load.py ../datasets/yeast_full.txt structures/yeast-AF2 saved_graphs
"""

import collections
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
from typing import Callable, List, Dict, Union
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
    verbose: bool = True,
) -> Dict[int, Dict[str, Union[str, sg.StellarGraph]]]:
    in_graphs: List[Dict[str, Union[str, nx.Graph]]] = list(graphs.values())

    out_graphs = {}

    for idx, graph in enumerate(in_graphs):
        g       = graph['graph'].copy()
        kinase  = graph['kinase']
        psite   = graph['psite']
        res     = graph['res']

        if g is None:
            pass

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
Get list of graph labels. 
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

    parent = os.path.dirname(in_path)
    save_path = os.path.join(parent, "embeddings")

    print(f"Input file is {in_path}.")

    # Unpickle
    infile = open(in_path,'rb')
    loaded_graphs = pickle.load(infile)
    infile.close()

    loaded = list(loaded_graphs.values())

    graphs = []
    graph_labels = []
    
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
        kinase = loaded[i]['kinase']
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

            #feature = [dist_to_psite, rsa, b_fac, *m]  # 10 long   
            feature = [dist_to_psite, rsa, b_fac, m[2], m[3], m[4], m[5]]  # 10 long 

            node_data["feature"] = feature
            
        # Create sg instance from graph, using `feature` vector. 
        g_attr = StellarGraph.from_networkx(g, node_features="feature")
        graphs.append(g_attr)
        graph_labels.append(kinase)

    
    #print(graphs[5].info())

    # Generator
    generator = sg.mapper.PaddedGraphGenerator(graphs)

    gc_model = sg.layer.GCNSupervisedGraphClassification(
        #[32, 16, 8], ["relu", "relu", "relu"], generator, pool_all_layers=True
        [32, 16], ["relu", "relu"], generator, pool_all_layers=True
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

    history = pair_model.fit(train_gen, epochs=epochs, verbose=0)
    #sg.utils.plot_history(history)

    embeddings = embedding_model.predict(generator.flow(graphs))

    print(f"Embeddings have shape {embeddings.shape}")

    # Save embeddings 
    print("Saving embeddings...")
    outfile =  open(save_path, 'wb')
    pickle.dump(embeddings, outfile)
    outfile.close()
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