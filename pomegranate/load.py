#! /usr/bin/env python3
"""Load graphs from list of PDB files"""


"""
usage: python load.py ../datasets/yeast_full.txt structures/yeast-AF2

"""


import pathlib
from graphein.protein.graphs import construct_graph
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.edges.distance import *	

from graphein.protein.config import DSSPConfig
from graphein.protein.features.nodes import rsa

from pathlib import Path

import click as c



import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'


import os
import urllib as ul

from definitions import STRUCTURE_PATH
from traitlets import default
from validate import get_database

# TODO: use Path types instead of strs
def load_graphs(
    pdb_path: str = None,       # directory containing local pdb files
    psite_list: str = None,          # path to file containing list of psites
    radius_threshold: float = 10.0,
    rsa_threshold: float = 0.8,
):

    psite_path = Path(psite_list)
    if not psite_path.is_file():
        raise ValueError(f"No such file {psite_list}")
    
    df = pd.read_csv(psite_path, sep='\t', header=0)

    
    accs = [a for a in df['acc']]

    # remove duplicates
    accs = list(set(accs))
    
    # Default directory
    if pdb_path == None:
        pdb_path = STRUCTURE_PATH + "/yeast"

    if not os.path.isdir(pdb_path):
        raise ValueError(f"Path {pdb_path} is not a directory.")
    
    # for each entry: 
    for acc in accs:

        filename = f"{pdb_path}/{acc}.pdb" 
        
        if not os.path.exists(filename):

            #print(f"No such file {filename}")
            print(f"Downloading {acc} from AF2...")

            filename = acc + ".pdb"
            url = f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v2.pdb"

            ul.request.urlretrieve(url, pdb_path+f"/{filename}")
    
    
        # Phosphosite 
        #subgraph = 
    pdb_dir = pdb_path 

    # Each phosphosite
    graphs = {}

    edge_fns = [
            add_aromatic_interactions,
            add_hydrophobic_interactions,
            add_aromatic_sulphur_interactions,
            add_cation_pi_interactions,
            add_disulfide_interactions,
            add_hydrogen_bond_interactions,
            add_ionic_interactions,
            add_peptide_bonds
        ]
    config = ProteinGraphConfig(
        edge_construction_functions=edge_fns, 
        graph_metadata_functions=[rsa], 
        dssp_config=DSSPConfig(),
        pdb_dir=pdb_dir,
    )
    for index, row in df.iterrows():

        #print(index, row['acc'], row['position'], row['code'])

        
        pdb_path = f"{pdb_dir}/{row['acc']}.pdb"
        g = construct_graph(config, pdb_path=pdb_path)      
        
        graphs[index] = g

        print(f"Graph {graphs[index].name} at {index}")

    return


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
@c.argument('phosphosite', nargs=1)
@c.argument('structures', nargs=1)
@c.option(
    "-r",
    "--radius",
    help="The threshold radius of the motif",
    type=float,
    default=10.0,
)
@c.option(
    "--rsa",
    "--rsa-threshold",
    help="The RSA threshold of the motif",
    type=float,
    default=0.5,
)
def main(
    phosphosite, 
    structures,
    radius,
    rsa,

    ):
    
    load_graphs(
        pdb_path = structures,
        psite_list = phosphosite,
        radius_threshold=radius,
        rsa_threshold=rsa,
    )

if __name__ == "__main__":
    main()