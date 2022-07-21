#! /usr/bin/env python3
"""Load graphs from list of PDB files"""


"""
usage: python load.py ../datasets/yeast_full.txt structures/yeast-AF2 saved_graphs
"""

#from distutils.log import debug
from cluster_graphs import nx_to_sg, GRAPH_NODE_FEATURES

import re

# Pomegranate
from protein.phosphosite import get_surface_motif
from validate import get_database

import pickle
import pathlib
from graphein.protein.graphs import construct_graph
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.edges.distance import *	

from graphein.utils.utils import protein_letters_3to1_all_caps as aa3to1

from graphein.protein.config import DSSPConfig
from graphein.protein.features.nodes import rsa

from pathlib import Path
from typing import List, Dict, Union, Tuple

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
    node_features: List[str] = None
    
) -> Dict[int, Dict[str, Union[str, nx.Graph]]]:    # TODO: use typed dict with each key/value pair being defined
    #graph = {'graph': g, 'kinase': kinase, 'psite': psite, 'res': res}
    # graphs[index] = graph

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

    # TODO 
    # Detect duplicates i.e. same kinase, same graph 
    # Option to store kinases as list instead of a whole new graph
    # Option to use 'known kinases' instead of 'unknown' -- as a CLI arg


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

    stats = dict(num_fail=0, num_success=0, num_skip=0)
    for idx, row in df.iterrows():

        #print(index, row['acc'], row['position'], row['code'])
        acc             = row['acc']
        res_pos         = row['position'] 
        res_code: str   = row['code']

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

                pos: int = int(res_pos)
                res: str = list(g.nodes())[pos-1]

                psite: Dict = g.nodes(data=True)[res]
                

                psite_res: str = str(psite['residue_name'])
                psite_num: int = int(psite['residue_number'])

                g = get_surface_motif(g, site=pos, r=radius_threshold, asa_threshold=rsa_threshold) 

                # Assert that phosphosite residue is same as what we expected 
                assert aa3to1(psite_res) == res_code, f"Residue mismatch {psite_res} and {res_code}"
                assert aa3to1(res.split(':')[1]) == res_code, f"Residue mismatch {res} and {res_code} {pos}"
                assert psite_num == pos

                g.name += f" @ {pos} {res_code}"


                
                # Assert that phosphosiste is included in the graph.  
                # TODO: display green on the terminal output if it is included; 
                # Display red on terminal if it is excluded (and --force was used.)

                graph = {'graph': g, 'kinase': kinase, 'psite': psite, 'res': res}
                graphs[index] = graph

                psite_contained = res in list(g.nodes())


                stats['num_success'] += 1
                
                if debug:
                    print(f"[{index}] Constructing graph from {acc}...", end=" ")
                    

                    if psite_contained: print('\x1b[6;30;42m' + '[PSITE]' + '\x1b[0m', end=" ")
                    else: print('\x1b[1;37;41m' + '[PSITE]' + '\x1b[0m', end=" ")
                    print(f"DONE.  Graph {graphs[index]['graph'].name} | PSITE: {res} | KINASE: {kinase}", end="")

                    #print(f"\t{'YES' if psite_contained else 'NO'}", end=" ")
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


"""
TODO: command line options

`-n` print what would be done i.e. dryrun
`-v` verbose 
`-q` quiet
`-d` debug 
`-i` interactive i.e. prompt
`-f` force, i.e. overwrite without asking etc.
`-c` config, i.e. provide file for describing how loading will work 


TODO: features

- autoskip (if ~enough PDB files are in directory, skip?)
- Warn user if skip used, but directory doesn't have any (many?) PDB files in it
- provide config file on how to load graphs 
- i.e. output file, edge features, node features (can do 'staging' i.e. before 
getting subgraph, we save temporarily so future config files can use the same 
saved graphs.  ), output format (e.g. nx, sg) 
"""

@c.command()
@c.argument('phosphosite', nargs=1)
@c.argument('structures', nargs=1)
@c.argument('graphs', nargs=1)
@c.option(
    "--verbose",
    "-v",
    is_flag=True,
    help="Show extensive program output."
)
@c.option(
    "--debug",
    "-d",
    is_flag=True,
    help="Show extensive program output for debugging."
)
@c.option(
    "--quiet",
    "-q",
    is_flag=True,
    help="Suppress program output."
)
@c.option(
    "--dry-run",
    "--dryrun",
    "-n",
    "is_dryrun",
    is_flag=True,
    help="Print out what the program would do without loading the graphs.",


)
@c.option(
    "--unique",
    "-u",
    is_flag=True,
    help="Only construct graphs for unique motifs.  Duplicate entries (i.e. with different kinases) are ignored."
)
@c.option(
    "--download/--skip-download", " /-S",
    help="Skip downloading protein structures into the supplied directory.  If the required structure does not exist, graph construction for that accession ID will be skipped.",
    default=True,
)
@c.option(
    # TODO: support multiple formats selected at once; i.e. saves more than one file.
    "-o",
    "--graph-format",
    "--output-format",
    "--format",
    type=c.Choice(['NetworkX', 'StellarGraph', 'nx', 'sg'], case_sensitive=False),
    help="Save graphs as NetworkX or StellarGraph instances with feature preprocessing. ",
    default="NetworkX",
    show_default=True,
)
@c.option(
    "-N",
    "--num-psites",
    help="Only consider the first N motifs in a dataset.  Graph construction will continue until N graphs are made, or the end of the dataset is reached.",
    type=c.INT,
    default=-1, 
)
@c.option(
    "-r",
    "--radius",
    help="The threshold radius of the motif",
    type=c.FLOAT,
    default=10.0,
    show_default=True,
)
@c.option(
    "--rsa",
    "--rsa-threshold",
    help="The RSA threshold of the motif",
    type=c.FLOAT,
    default=0.0,
    show_default=True,
)
@c.option(
    "--node-features",
    "--nf",
    is_flag=False,
    default=','.join(GRAPH_NODE_FEATURES), show_default=True,
    metavar="<node_features>",
    type=c.STRING,
    help="Which node features to include in the constructed graphs."
)
@c.option(
    "--edge-features",
    "--ef",
    is_flag=False,
)
@c.option(
    "--config",
    "-c",
    help="Path to config.yml file used to specify how to construct the graphs.",
    # TODO: check if path right here?
    default="config.yml",
    show_default=True,
)
def main(
    # POSITIONAL ARGUMENTS:
    phosphosite,    
    structures,
    graphs,

    # FLAGS:
    verbose,
    debug,
    quiet,
    is_dryrun,
    unique,
    download,

    # PARAMETERS:
    graph_format,
    num_psites,
    radius,
    rsa,
    node_features,
    edge_features,
    config,
):
    # TODO: check that the phosphosite used in graphein construction IS IN FACT
    # the same residue that we see from the phosphosite.  Check for 'off by 1' error?


    # TODO: check for duplicates (i.e. same kinase, same graph motif position )

    verbose = True # TODO: Remove this


    if debug: 
        verbose = True

    if is_dryrun:
        verbose = False

    

    
    
    if graph_format.lower() in ["sg", "stellargraph"]:
        output_format = "sg"
    elif graph_format.lower() in ["networkx", "nx"]:
        output_format = "nx"
    else:
        raise NotImplementedError(f"Output format '{output_format}' not implemented.")

    
    # TODO: ensure that psite is always included; regardless of RSA
    # TODO: might be multiple identical graphs (i.e. only different thing in entry is the kinase)
    # have option to only count unique graphs.  Or utilise different entries. 

    # TODO: 'stage' the preprocessing of graphs by first. 
    # 1.  load set of graphs with given n features (edges etc.)
    # 2. 'select' subset of features from already saved list of graphs
    # 3.  i.e. have a saved version of "all graphs with edge features [a,b,c]" 
    #           saved version of "all graphs with edge features [d, e]" etc.
    # 4. load from this; then save all graphs with radius 10A, 0.2 RSA or whatever.  
    # 5.  so entire loading doesn't have to occur again; when a different radius is chosen. 

    # 6.  split clustering into visualisation and clustering; so we don't have to perform GCN 
    # if we just want to try a diff way of viewing with tSNE or something. 

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

    print("\nGRAPH LOADER\n------------")
    
    if verbose: print(f"Output file is {out_path}.")

    # Handle features
    node_features = [f.strip() for f in node_features.split(',')]
    #edge_features = [f.strip() for f in edge_features.split(',')]

    if verbose: print(f"Using node features: {node_features}")

    
    if is_dryrun:
        n = num_psites if num_psites > 0 else "all"
        download_status = f"Will attempt to download PDB files into {structures}" if download else "Skipping download from AF2"
        print(f"""
        {download_status}

        Output file is {out_path}.

        Create {n} graphs from dataset {phosphosite}. 
        Radius threshold:\t{radius}
        RSA threshold:\t{rsa}
        Format:\t{graph_format}

        Features:
        \tNodes:\t{node_features}
        \tEdges:\t


    
        """)

        raise NotImplementedError(f"--dryrun not implemented yet. ")
        exit
    
    

    graphs = load_graphs(
        pdb_path = structures,
        psite_list = phosphosite,
        radius_threshold=radius,
        rsa_threshold=rsa,
        num_psites=num_psites,
        download=download,
        node_features=node_features
    )

    print(f"Created {len(graphs.values())} graphs with radius {radius} and RSA {rsa}")


    # Save data
    data = {}

    if output_format == "sg":
        
        
        if verbose: print(f"Converting graphs to StellarGraph instances...", end=" ") 

        graphs = nx_to_sg(graphs, include_features=node_features)

        if verbose: print(f"DONE.")

        if debug:
            g = graphs[0]['graph']
            print(f"Graph 0 is type {type(g)}")
            print(g.info())

    elif output_format == "nx":
        if verbose: print(f"Graphs are in NetworkX format.")
    else:
        raise NotImplementedError(f"Output format '{output_format}' not implemented.")
    
    data['graphs_dict'] = graphs
    data['format'] = output_format

    # Save graphs to file
    if verbose: print(f"Saving graphs to {out_path} ...", end=" ")
    outfile =  open(out_path, 'wb')
    pickle.dump(data, outfile)
    outfile.close()
    if verbose: print("DONE.")



    return

 
if __name__ == "__main__":
    main()
