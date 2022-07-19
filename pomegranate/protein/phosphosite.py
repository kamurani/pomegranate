### EXPORTS

__all__ = ['get_protein_graph', 'get_protein_subgraph_radius', 
            'get_adjacency_matrix_plot',
            'get_phosphosites']




### DEFINITIONS
from definitions import ROOT_DIR, STRUCTURE_PATH

### External libraries
import os
import numpy as np

# Plotly 
import plotly.express as px

from graphein.protein.config import ProteinGraphConfig
from graphein.protein.graphs import construct_graph
from graphein.protein.utils import download_alphafold_structure

from graphein.protein.visualisation import plot_distance_matrix, plotly_protein_structure_graph, plot_chord_diagram	
from graphein.protein.subgraphs import extract_subgraph_from_point, extract_k_hop_subgraph, extract_subgraph, extract_surface_subgraph


from graphein.protein.edges.distance import *	

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
def get_phosphosites(g, residues=['SER', 'THR', 'TYR', 'HIS']):
    
    psites = extract_subgraph(g, 
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
            add_aromatic_interactions,
            add_hydrophobic_interactions,
            add_aromatic_sulphur_interactions,
            add_cation_pi_interactions,
            add_disulfide_interactions,
            add_hydrogen_bond_interactions,
            add_ionic_interactions,
            add_peptide_bonds
            ]

        # Use structure path of already downloaded PDB file (if it exists) for DSSP calculation.
        pdb_path = STRUCTURE_PATH + '/' + id + '.pdb'

        from graphein.protein.config import DSSPConfig
        from graphein.protein.features.nodes import rsa
        config = ProteinGraphConfig(edge_construction_functions=edge_fns, 
                                    graph_metadata_functions=[rsa], 
                                    dssp_config=DSSPConfig(),
                                    pdb_path=pdb_path,
        )   
    

    protein_path = id + '.pdb'

    if database in ['AlphaFold', 'SWISS_PROT']:
        #NOTE: Might have to remove SWISS_PROT. Not all SP have AF structures
        protein_path = download_alphafold_structure(id, aligned_score=False, out_dir=STRUCTURE_PATH)

    # if use_alphafold:
    #     pdb_path = download_alphafold_structure(id, aligned_score=False, out_dir=STRUCTURE_PATH)
   
    # TODO: separate structures into alphafold / pdb. 
    # Check if this file has been downloaded before.
    # if os.path.isfile(pdb_path):
    #     print(f"Using local PDB file for {id}.")
    #     g = construct_graph(config=config, pdb_path=pdb_path)
    if os.path.isfile(protein_path):
        print(f"Using local file for {id}.")
        g = construct_graph(config=config, pdb_path=protein_path)
    else:
        print(f"Retrieving {id}...")
        g = construct_graph(config=config, pdb_code=id)

    
    # TODO: check if file exists and download if not. 
   
    
    # construct graph
   
    return g

'''
Given graph ``g`` get subgraph within radius of psite, and surface residues above 
ASA threshold. 
'''
def get_surface_motif(g=None, site=1, r=10, asa_threshold=0.5):

    
    s_g = get_protein_subgraph_radius(g=g, site=site, r=r)

    if asa_threshold:
        try:
            surface = extract_surface_subgraph(s_g, asa_threshold, 
                                                recompute_distmat=True,
                                                filter_dataframe=True
            )
        except:
            raise ValueError("Specified graph does not have RSA metadata.")
        return surface
    else:
        return s_g # Don't consider surface if asa is None

    
    
'''
Given a graph ``g`` get a subgraph from radius and known phos site
'''
def get_protein_subgraph_radius(g=None, site=1, r=10):
   
    if isinstance(site, str):
        try:
            site = int(site.split(':')[-1])
        except ValueError:
            raise ValueError("Specified phospho site isn't in correct format.")            
        
    # get centre point
    #index = query['phosphosite'] - 1
    index = site - 1
    phos_point = tuple(g.graph['coords'][index])
    
    # Get subgraph
    s_g = extract_subgraph_from_point(g, centre_point=phos_point, 
                                        radius=r, 
                                        recompute_distmat=True, 
                                        filter_dataframe=True)
    
    return s_g
    

# TODO: make separate function for extracting subgraph that is independent of the get_plot


def get_adjacency_matrix_plot(g=None, psite=1, title=None, order='seq'):
    
    # subgraph
    
    #fig = plot_distance_matrix(g, title=title)
    fig = motif_plot_distance_matrix(g, psite=psite, title=title, aa_order=order)
    return fig
