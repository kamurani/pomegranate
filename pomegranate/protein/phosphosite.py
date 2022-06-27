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
from graphein.protein.subgraphs import extract_subgraph_from_point, extract_k_hop_subgraph, extract_subgraph
	
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
def get_protein_graph(id=None, path=None, database="AlphaFold", out_dir='/tmp'):
    

    # TODO: option that 'saves' the PDB file instead of downloading it to /tmp 

    database = database.lower()

    config = ProteinGraphConfig()
    #query = dict(prot_id="Q9Y2X7", phosphosite=224)


    # protein id given; check for local copy
    if id:
        protein_path = STRUCTURE_PATH + '/' + id + '.pdb'
    # user-specified local file
    if path: 
        protein_path = path
        if not os.path.isfile(protein_path):
            raise ValueError("Specified path is not a PDB file.")
 
    if protein_path is None:
        raise ValueError("Must specify a protein by providing a PDB file path or protein code.")
    
    # TODO: separate structures into alphafold / pdb. 


    # Check if this file has been downloaded before.
    if protein_path and os.path.isfile(protein_path):
        print(f"Using local PDB file for {id}.")
        g = construct_graph(config=config, pdb_path=protein_path, pdb_code=id)
    else:
        if database == 'alphafold':
            protein_path = download_alphafold_structure(id, aligned_score=False, out_dir=out_dir)
            print(f"Retrieving {id} from AlphaFold...")
            g = construct_graph(config=config, pdb_path=protein_path)

        elif database == 'pdb':   
            try:
                g = construct_graph(config=config, pdb_code=id)
            except:
                raise ValueError(f"Invalid PDB code.")
            print(f"Retrieving {id} from PDB...")
    

    return g
    
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


def get_adjacency_matrix_plot(g=None, psite=1, title=None):
    
    # subgraph
    
    #fig = plot_distance_matrix(g, title=title)
    fig = motif_plot_distance_matrix(g, psite=psite, title=title)
    return fig
