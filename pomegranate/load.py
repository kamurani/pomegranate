#! /usr/bin/env python3
"""Load graphs from list of PDB files"""


"""
usage: python load.py ../datasets/yeast_full.txt structures/yeast-AF2 saved_graphs
"""

import re

# Pomegranate
from protein.phosphosite import get_surface_motif
from validate import get_database

import pickle
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
    rsa_threshold: float = 0.0,
    download: bool = True,
    num_psites: int = -1, # Make 1000 graphs as default
    verbose: bool = True,
    debug: bool = True,
    
):

    if debug:
        verbose = True


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
    
    if not download:
        if verbose:
            print(f"Skipping downloading.  Using {pdb_path} as structure directory...", end=" ")

        filenames = os.listdir(pdb_path)
        num_matches = 0 
        for f in filenames:
            if re.search(".pdb$", f):
                num_matches += 1
        if verbose:
            print(f"Found {num_matches} PDB files.", end=" ")
        
        if num_matches == 0:
            if verbose:
                print("Exiting...", end=" ")
            exit
        if verbose:
            print()
        
        
    else:
        # for each entry: 
        for acc in accs:

            filename = f"{pdb_path}/{acc}.pdb" 
            if os.path.exists(filename):
                if verbose:
                    print(f"{filename} already exists.")

            else:

                #print(f"No such file {filename}")
                

                filename = acc + ".pdb"

                if get_database(acc) == 'uniprot':
                    print(f"Downloading {acc} from AF2...", end=" ")
                    try:
                        url = f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v2.pdb"
                        ul.request.urlretrieve(url, pdb_path+f"/{filename}")
                        if verbose:
                            print("DOWNLOADED.")
                    except:
                        if verbose:
                            print(f"FAILED to download from AF2.")
                else:
                    if verbose:
                        print(f"Skipping non-uniprot ID '{acc}'.")
        
        
        # Phosphosite 
        #subgraph = 
    pdb_dir = pdb_path 
    # Each phosphosite
    graphs = {}

    edge_fns = [
            #add_aromatic_interactions,
            add_hydrophobic_interactions,
            #add_aromatic_sulphur_interactions,
            #add_cation_pi_interactions,
            #add_disulfide_interactions,
            #add_hydrogen_bond_interactions,
            #add_ionic_interactions,
            add_peptide_bonds
        ]
    config = ProteinGraphConfig(
        edge_construction_functions=edge_fns, 
        graph_metadata_functions=[rsa], 
        dssp_config=DSSPConfig(),
        pdb_dir=pdb_dir,
    )

    stats = dict(num_fail=0, num_success=0, num_skip=0)
    for idx, row in df.iterrows():

        #print(index, row['acc'], row['position'], row['code'])
        acc         = row['acc']
        pos         = row['position'] 
        res_code    = row['code']

        if type(row['kinases']) == str:
            kinase = row['kinases']
        else:
            kinase = "UNKNOWN"
        #print(f"row kinase: {type(row['kinases'])}")
        index = idx
        #if True:
        
        
        pdb_path = f"{pdb_dir}/{acc}.pdb"
        path = Path(pdb_path)
        if not path.is_file():
            stats['num_skip'] += 1
            if verbose:
                print(f"[{index}] No structure file {pdb_path} -- Skipping...")
        else:
            if verbose:
                print(f"[{index}] Constructing graph from {acc}...", end=" ")
            
            try:
                g = construct_graph(config, pdb_path=pdb_path) 
                res = list(g.nodes())[pos-1]

                psite = g.nodes(data=True)[res]
                g = get_surface_motif(g, site=pos) # use default thresholds
                g.name += f" @ {pos} {res_code}"

                graph = {'graph': g, 'kinase': kinase, 'psite': psite, 'res': res}
                graphs[index] = graph
                
                stats['num_success'] += 1
                
                if debug:
                    print(f"[{index}] Constructing graph from {acc}...", end=" ")
                    print(f"DONE.  Graph {graphs[index]['graph'].name}| PSITE: {res} | KINASE: {kinase}", end="")
                if verbose:
                    print("")
                
            except:
                graphs[index] = None
                stats['num_fail'] += 1
                if verbose:
                    print(f"FAILED.")
        
        # Exit if we have reached N graphs, and the user supplied an N
        if stats['num_success'] >= num_psites and num_psites > 0:
            break
            
    # Print stats      
    if verbose:
        print("STATISTICS")
        print("----------")
        print(f"{idx+1} phosphosites examined")
        print(f"{stats['num_skip']} graphs skipped")
        print(f"{stats['num_success']} graphs successfully constructed")
        print(f"{stats['num_fail']} graph constructions failed")
        print("")

    return graphs


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
@c.argument('graphs', nargs=1)
@c.option(
    "--download/--skip-download",
    help="Skip downloading protein structures into the supplied directory.  If the required structure does not exist, graph construction for that accession ID will be skipped.",
    default=True,
)
@c.option(
    "-N",
    "--num-psites",
    help="Only consider the first N motifs in a dataset.  Graph construction will continue until N graphs are made, or the end of the dataset is reached.",
    type=int,
    default=-1, 
)
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
    default=0.0,
)
def main(
    phosphosite, 
    structures,
    graphs,
    download,
    num_psites,
    radius,
    rsa,
    ):

    # TODO: ensure that psite is always included; regardless of RSA
    # TODO: might be multiple identical graphs (i.e. only different thing in entry is the kinase)
    # have option to only count unique graphs.  Or utilise different entries. 

    graph_path = Path(graphs)

    if graph_path.is_file():
        out_path = graph_path

        # TODO: check for overwrite if file exists.  Prompt. 
        # '--force, -f' flag to overwrite. 
    elif graph_path.is_dir():
        filename = "graph_objects"
        out_path = os.path.join(graph_path, filename)
    else:
        raise ValueError(f"No such directory {graph_path}")
    
     
    # TODO: check if filename exists.  Prompt for new one / overwrite. 

    print(f"Output file is {out_path}.")
    

    graphs = load_graphs(
        pdb_path = structures,
        psite_list = phosphosite,
        radius_threshold=radius,
        rsa_threshold=rsa,
        num_psites=num_psites,
        download=download,
    )

    print(f"Created {len(graphs.values())} graphs with radius {radius} and RSA {rsa}")

    # Save graphs to file
    print(f"Saving graphs to {out_path} ...", end=" ")
    outfile =  open(out_path, 'wb')
    pickle.dump(graphs, outfile)
    outfile.close()
    print("DONE.")



    return

 
if __name__ == "__main__":
    main()