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
@c.argument('graphs', nargs=1)
#@c.argument('structures', nargs=1)
#@c.argument('graphs', nargs=1)
def main(
    graphs,
):

    graph_path = Path(graphs)
    if not graph_path.is_dir():
        raise ValueError(f"No such directory {graph_path}")
    
     
    filename = "graph_objects"
    in_path = os.path.join(graph_path, filename)

    print(f"Input file is {in_path}.")

    # Unpickle
    infile = open(in_path,'rb')
    loaded_graphs = pickle.load(infile)
    infile.close()

    print(loaded_graphs[1])


if __name__ == "__main__":
    main()