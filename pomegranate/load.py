#! /usr/bin/env python3
"""Load graphs from list of PDB files"""


import pathlib
from graphein.protein.graphs import construct_graph
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.edges.distance import *	

from pathlib import Path

import click as c



import pandas as pd

import os
import urllib as ul

from definitions import STRUCTURE_PATH
from validate import get_database

# TODO: use Path types instead of strs
def load_graphs(
    pdb_path: str = None,       # directory containing local pdb files
    psite_list: str = None,          # path to file containing list of psites
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
def main(
    phosphosite, 
    structures,
    radius,
    rsa,
    
    ):
    
    load_graphs(
        pdb_path = structures,
        psite_list = phosphosite,
    )

if __name__ == "__main__":
    main()