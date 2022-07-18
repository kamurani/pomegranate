"""Web application for testing."""

# Run this app with `python app_testing.py` and
# visit http://127.0.0.1:8051/ in your web browser.
from visualisation.plot import motif_asteroid_plot
from definitions import STRUCTURE_PATH

import dash
from dash import Dash, dcc, html, Input, Output

from graphein.protein.graphs import construct_graph
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.edges.distance import *	

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = Dash(__name__, external_stylesheets=external_stylesheets)

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
#pdb_path = STRUCTURE_PATH + '/' + id + '.pdb'

pdb_path = STRUCTURE_PATH+"/Q9Y2X7.pdb"

from graphein.protein.config import DSSPConfig
from graphein.protein.features.nodes import rsa
config = ProteinGraphConfig(edge_construction_functions=edge_fns, 
                            graph_metadata_functions=[rsa], 
                            dssp_config=DSSPConfig(),
                        
                            pdb_dir=STRUCTURE_PATH,
)
g = construct_graph(config, pdb_path=pdb_path)

#g = construct_graph(config, pdb_path="/home/cam/DESN2000/pomegranate/structures"+"/Q9Y2X7.pdb")


psite = list(g.nodes)[200]
#psite = 'A:SER:498'
fig = motif_asteroid_plot(
    g=g, 
    node_id=psite,
    size_nodes_by="rsa",
    node_size_multiplier=80,
    colour_nodes_by="hydrophobicity",
)

'''
fig2 = motif_asteroid_plot(
    g=g, 
    node_id=psite,
    size_nodes_by="degree",
)
'''

app.layout = html.Div([
    dcc.Graph(figure=fig),
    #dcc.Graph(figure=fig2),
])



def run():
    app.run_server(debug=True, port=8051)

if __name__ == '__main__':
    run()

