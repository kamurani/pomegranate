### EXPORTS

__all__ = ['get_protein_graph', 'get_protein_subgraph_radius', 
            'get_adjacency_matrix_plot',
            'get_phosphosites']



# Thr12 on Ubiquitin is phosphorylated (5nvg.pdb)

from graphein.protein.config import ProteinGraphConfig
from graphein.protein.graphs import construct_graph
from graphein.protein.utils import download_alphafold_structure

from graphein.protein.visualisation import plot_distance_matrix, plotly_protein_structure_graph, plot_chord_diagram	
from graphein.protein.subgraphs import extract_subgraph_from_point, extract_k_hop_subgraph, extract_subgraph
	
# Custom plot function
from visualisation.plot import motif_plot_distance_matrix

import numpy as np



'''
TODO: 
- colour by similarity (AA) instead of distance
'''



# Plotly 
import plotly.express as px


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
def get_protein_graph(id=None):
    
    config = ProteinGraphConfig()
    #query = dict(prot_id="Q9Y2X7", phosphosite=224)
    
    out_dir = './structures' 
    
    protein_path = out_dir + '/' + id + '.pdb'
    
    # TODO: check if file exists and download if not. 
    #protein_path = download_alphafold_structure(prot_id, aligned_score=False, out_dir=out_dir)
    
    # construct graph
    g = construct_graph(config=config, pdb_path=protein_path)
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
