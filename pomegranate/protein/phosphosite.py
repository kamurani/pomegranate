### EXPORTS

__all__ = ['get_protein_graph', 'get_protein_subgraph_radius', 
            'get_adjacency_matrix_plot',
            'get_phosphosites']




### DEFINITIONS
from definitions import ROOT_DIR, STRUCTURE_PATH, SAVED_PDB_DIR, SAVED_GRAPHS_DIR, AF_VERSION

### External libraries
import os
import numpy as np
import json

# Plotly 
import plotly.express as px

from graphein.protein.config import ProteinGraphConfig
from graphein.protein.graphs import construct_graph
from graphein.protein.utils import download_alphafold_structure

from graphein.protein.visualisation import plot_distance_matrix, plotly_protein_structure_graph, plot_chord_diagram	
from graphein.protein.subgraphs import extract_subgraph_from_point, extract_k_hop_subgraph, extract_subgraph, extract_surface_subgraph


import graphein.protein.edges.distance as g_dist

# Custom plot function
from visualisation.plot import motif_plot_distance_matrix

'''
TODO: 
- colour by similarity (AA) instead of distance
'''

# TODO convert from single to triple letter codes

'''
Return phosphosites (sorted)
'''
def get_phosphosites(g, residues=['SER', 'THR', 'TYR', 'HIS'], rsa_threshold=0.5):
    
    surface_nodes = extract_subgraph(g, rsa_threshold=rsa_threshold)
    psites = extract_subgraph(surface_nodes,
                            return_node_list=True,
                            residue_types=residues)
    
    psites_sorted = sorted(psites, key=lambda x: int(x.split(':')[-1]))
    return psites_sorted
    
# TODO: make this function receive a `list` of dict(id=id, site=site) 
# this function then returns a list of graphs
def get_protein_graph(id=None, config=None, database='PDB'):
    # Graph configuration
    if not config:
        config = ProteinGraphConfig()   # default graph config file from graphein
    
    if config in ["asa", "rsa"]:

        # Edge functions
        edge_fns = [
            g_dist.add_aromatic_interactions,
            g_dist.add_hydrophobic_interactions,
            g_dist.add_aromatic_sulphur_interactions,
            g_dist.add_cation_pi_interactions,
            g_dist.add_disulfide_interactions,
            g_dist.add_hydrogen_bond_interactions,
            g_dist.add_ionic_interactions,
            g_dist.add_peptide_bonds
            ]

        # Use structure path of already downloaded PDB file (if it exists) for DSSP calculation
        pdb_path = STRUCTURE_PATH + '/' + id + '.pdb'

        from graphein.protein.config import DSSPConfig
        from graphein.protein.features.nodes import rsa
        config = ProteinGraphConfig(edge_construction_functions=edge_fns, 
                                    graph_metadata_functions=[rsa], 
                                    dssp_config=DSSPConfig(),
                                    pdb_path=pdb_path,
        )   
    
    # NOTE: File paths use '\' in windows systems
    # NOTE: Need different prot_dir for each DB
    protein_path = SAVED_PDB_DIR + id + '.pdb'
    if database in ['AlphaFold', 'SWISS_PROT']:
        #NOTE: Might have to remove SWISS_PROT. Not all SP have AF structures
        print("AF or SP")
        # NOTE: DEBUGGING: construct graph is always looking at examples/pdbs. Just put everything in examples/pdbs
        # protein_path = download_alphafold_structure(id, aligned_score=False, out_dir=STRUCTURE_PATH)
        protein_path = download_alphafold_structure(id, AF_VERSION, aligned_score=False, out_dir=SAVED_PDB_DIR)
        #protein_path = STRUCTURE_PATH + id + '.pdb'
        print("After")

    # if use_alphafold:
    #     pdb_path = download_alphafold_structure(id, aligned_score=False, out_dir=STRUCTURE_PATH)
   
    # TODO: separate structures into alphafold / pdb. 
    # Check if this file has been downloaded before.
    # if os.path.isfile(pdb_path):
    #     print(f"Using local PDB file for {id}.")
    #     g = construct_graph(config=config, pdb_path=pdb_path)

    # Check if graph exists
    graph_path = f'{SAVED_GRAPHS_DIR}/{id}_{database}.json'
    if os.path.isfile(graph_path):
        with open(graph_path, "r") as f:
            print(f"Using local graph for {id} from {database}")
            g = json.load(f)
    else:
        # Graph doesn't exist
        if os.path.isfile(protein_path):
            print(f"Using local file for {id}.")
            #g = construct_graph(config=config, pdb_path=protein_path)
            if database == 'PDB':
                g = construct_graph(config=config, pdb_path=protein_path)
            else:
                g = construct_graph(config=config, pdb_path=protein_path)
        else:
            print(f"Retrieving {id}...")
            if database == 'PDB':
                #g = construct_graph(config=config, pdb_code=id)
                g = construct_graph(config=config, pdb_code=id)
            else: # NOTE: FIX THIS. BAD STYLE. Same line as 119
                #g = construct_graph(config=config, pdb_path=protein_path)
                g = construct_graph(config=config, pdb_path=protein_path)

    # TODO: check if file exists and download if not. 
   
    
    # construct graph
   
    return g

'''
Given graph ``g`` get subgraph within radius of psite, and surface residues above 
ASA threshold. 
'''
def get_surface_motif(g, site, r=10, asa_threshold=0.5):

    s_g = get_protein_subgraph_radius(g=g, site=site, r=r)

    if asa_threshold:
        try:
            surface = extract_surface_subgraph(s_g,
                                               asa_threshold,
                                               recompute_distmat=True,
                                               filter_dataframe=True)
        except:
            raise ValueError("Specified graph does not have RSA metadata.")

        #surface.add_node(psite_node) # Restore psite node if it was removed
        return surface
    else:
        #s_g.nodes(data=True)[res] = psite_node
        return s_g # Don't consider surface if asa is None

    
    
'''
Given a graph ``g`` get a subgraph from radius and known phos site
'''
def get_protein_subgraph_radius(g, site, r=10.0):

    # get centre point   
    try:
        x_y_z = node_coords(g, site)
    except ValueError:
        raise ValueError("Specified phospho site isn't in correct format.")       
    
    # Get subgraph
    s_g = extract_subgraph_from_point(g, centre_point=x_y_z, radius=r)
    
    return s_g
    

# TODO: make separate function for extracting subgraph that is independent of the get_plot


def get_adjacency_matrix_plot(g=None, psite=1, title=None, order='seq', colour='viridis_r'):
    
    # subgraph
    
    #fig = plot_distance_matrix(g, title=title)
    fig = motif_plot_distance_matrix(g, psite=psite, title=title, aa_order=order, colour=colour)
    return fig

# Get x, y, z coordinates from a given node in a graph
# Input: - Graph g
#        - str node (e.g. 'A:ARG:1')
# Output: tuple (x, y, z)
def node_coords(g, node):

    df = g.graph['pdb_df']

    coords = df.loc[df.node_id == node][['x_coord','y_coord','z_coord']]
    if not coords.empty:
        return coords.values[0]
    else:
        return None