import os
from dash import Dash, dcc, html, callback, Input, Output, dependencies
import plotly.express as px
import pandas as pd
import numpy as np
import math

from protein.phosphosite import get_surface_motif
from protein.interactions import add_distance_threshold
from visualisation.plot import motif_asteroid_plot
from json_conversion.graphein_to_json import g_to_json, load_prot_graph

import os
from xml.etree.ElementInclude import include
from dash import Dash, html, dcc, Input, Output
import pandas as pd
import plotly.express as px


from definitions import EMBEDDINGS_PATH, STRUCTURE_HUMAN_PATH, SAVED_CLUSTER_GRAPHS_PATH
from utils.amino_acid import aa1letter
from visualisation.plot import motif_plot_distance_matrix
from protein.phosphosite import get_protein_graph

from graphein.protein.visualisation import plot_distance_matrix
from graphein.utils.utils import protein_letters_3to1_all_caps as aa3to1

from graphein.protein.config import ProteinGraphConfig
from graphein.protein.graphs import construct_graph




import networkx as nx


from typing import Callable, Dict, List, Union


"""
Load from save
"""
def load_graphs(

):
    pass



"""
Construct graphs and serialize
"""
def construct_graphs(
    df,
    pdb_dir: str, 
    out_dir: str = SAVED_CLUSTER_GRAPHS_PATH,
    overwrite: bool = False, 

):

    """
    :param overwrite: If ``True``, ignore already existing files in ``out_dir`` with the same name and overwrite.
    :type overwrite: bool
    
    """

    if not os.path.isdir(pdb_dir):
        raise ValueError(f"No specified directory for structures '{pdb_dir}'")

    if not os.path.isdir(out_dir):
        raise ValueError(f"No specified directory for JSON graphs '{out_dir}'")

    print(f"Saving to {out_dir}")

    dff = df
    dff["Protein ID"] = dff["Protein ID"].apply(    # get just protein name
        lambda x: pd.Series(str(x).split('@')[0])
    )
    dff["Psite ID"] = df["Protein ID"] + df["Phosphosite"]  # concatenate columns
                             # unique protein / psite pairs

    assert os.path.isdir(pdb_dir)

    graphs: Dict[str, nx.Graph] = {}
    
    config = ProteinGraphConfig()

    database = "AF2" # TODO: in future, this may be in the `df` to check for. 

    for prot, psite in [p.split() for p in dff['Psite ID'].unique()]:
        
        name = f"{prot}_{psite}"

        pdb_path = f"{pdb_dir}/{prot}.pdb"

        filename = name + f"_{database}.json"
        
        path = os.path.join(out_dir, filename)
        #print(f"path: {path}")
        if not overwrite and os.path.isfile(path):
            print(f"Skipping {name} ...")
            continue

        #print(f"Constructing {name} ...")

        try: 
            g = construct_graph(config, pdb_path=pdb_path) 

            """
            Note: ideally we would be obtaining whole graph to be used later; but for performance 
            we are restricting clustering tab visualisations of adj matrices etc. to use radius <= 10A 

            Further subgraphs can be selected dynamically by other functions.
            """
            g = get_surface_motif(g, site=psite, r=10, asa_threshold=0) 

        except: 
            g = {}


        g.name = name

        g_to_json(
            g, 
            prot_id=name,
            db_name="AF2",
            save_path=out_dir,
        )

        print(f"Saved.")
        
        
        #graphs[name] = g    # DON'T store in memory for now; save to json file


        # TODO: only store adjacency matrix? 

    return graphs


