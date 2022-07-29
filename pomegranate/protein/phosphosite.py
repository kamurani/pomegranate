### EXPORTS

__all__ = ['get_protein_graph', 'get_protein_subgraph_radius', 
            'get_adjacency_matrix_plot',
            'get_phosphosites']




### DEFINITIONS
from definitions import ROOT_DIR, STRUCTURE_PATH

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
    prot_dir = '../examples/pdbs/'
    protein_path = prot_dir + id + '.pdb'

    if database in ['AlphaFold', 'SWISS_PROT']:
        #NOTE: Might have to remove SWISS_PROT. Not all SP have AF structures
        print("AF or SP")
        protein_path = download_alphafold_structure(id, aligned_score=False, out_dir=STRUCTURE_PATH)
        print("After")

    # if use_alphafold:
    #     pdb_path = download_alphafold_structure(id, aligned_score=False, out_dir=STRUCTURE_PATH)
   
    # TODO: separate structures into alphafold / pdb. 
    # Check if this file has been downloaded before.
    # if os.path.isfile(pdb_path):
    #     print(f"Using local PDB file for {id}.")
    #     g = construct_graph(config=config, pdb_path=pdb_path)

    # Check if graph exists
    graph_dir= '../graphs'
    graph_path = f'{graph_dir}/{id}_{database}.json'
    if os.path.isfile(graph_path):
        with open(graph_path, "r") as f:
            print(f"Using local graph for {id} from {database}")
            g = json.load(f)
    else:
        # Graph doesn't exist
        if os.path.isfile(protein_path):
            print(f"Using local file for {id}.")
            g = construct_graph(config=config, pdb_path=protein_path)
        else:
            print(f"Retrieving {id}...")
            if database == 'PDB':
                g = construct_graph(config=config, pdb_code=id)
            else: # NOTE: FIX THIS. BAD STYLE. Same line as 119
                g = construct_graph(config=config, pdb_path=protein_path)

    # TODO: check if file exists and download if not. 
   
    
    # construct graph
   
    return g

'''
Given graph ``g`` get subgraph within radius of psite, and surface residues above 
ASA threshold. 
'''
def get_surface_motif(
    g: nx.Graph = None, 
    site: Union[int, str] = 1, 
    r: float = 10.0, 
    asa_threshold: float = 0.5,
):
    # res = list(g.nodes())[site-1]
    # # print("res is", res)
    # psite_node = g.nodes(data=True)[res]
    
    s_g = get_protein_subgraph_radius(g=g, site=site, r=r)



    if asa_threshold:
        try:
            surface = extract_surface_subgraph(s_g, asa_threshold, 
                                                recompute_distmat=True,
                                                filter_dataframe=True
            )
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
def get_protein_subgraph_radius(
    g: nx.Graph = None, 
    site: Union[int, str] = 1, 
    r: float = 10.0,
):
   
    if isinstance(site, str):
        try:
            site = int(site.split(':')[-1])
        except ValueError:
            raise ValueError("Specified phospho site isn't in correct format.")            
        
    # get centre point
    #index = query['phosphosite'] - 1
    index = site - 1
    phos_point = np.array(g.graph['coords'][index]) #TODO: TURN BACK TO TUPLE IF LIST DOEESN'T WORK

    print(phos_point)

    
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
