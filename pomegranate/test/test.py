from graphein.protein.graphs import construct_graph
from graphein.protein.utils import download_alphafold_structure
import graphein.protein.edges.distance as g_dist
from graphein.protein.config import ProteinGraphConfig
import networkx.readwrite as nx
import json

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
pdb_path = ''

from graphein.protein.config import DSSPConfig
from graphein.protein.features.nodes import rsa
config = ProteinGraphConfig(edge_construction_functions=edge_fns, 
                            graph_metadata_functions=[rsa], 
                            dssp_config=DSSPConfig(),
                            pdb_path=pdb_path,
)   
id = "q5vsl9"
protein_path = download_alphafold_structure(id, aligned_score=False, out_dir='')
g = construct_graph(pdb_path=protein_path)
del g.graph["config"]

for k, v in g.graph.items():
    try:
        g.graph[k] = v.to_json()
    except AttributeError:
        try:
            g.graph[k] = v.tolist()
        except AttributeError:
            continue

# Ensure node data is in JSON format
for n, d in g.nodes(data=True):
    for k, v in d.items():
        try:
            d[k] = v.to_json()
        except AttributeError:
            try:
                d[k] = v.tolist()
            except AttributeError:
                continue

# Ensure edge data is in JSON format
for _, _, d in g.edges(data=True):
    for k, v in d.items():
        try:
            d[k] = v.to_json()
        except AttributeError:
            try:
                d[k] = v.tolist()
            except AttributeError:
                try:
                    d[k] = list(v)
                except AttributeError:
                    continue

with open("test.json", 'w') as f:
    # tmp = json_graph.node_link_data(g)
    # f.write(json.dumps(tmp))
    tmp = nx.json_graph.node_link_data(g)
    f.write(json.dumps(tmp))