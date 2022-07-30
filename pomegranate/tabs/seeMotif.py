
from dash import Dash, dcc, html, callback, Input, Output, dependencies
import plotly.express as px
import pandas as pd
import numpy as np
import math

import networkx.readwrite as nx
import json

from protein.phosphosite import *
from protein.phosphosite import get_surface_motif
from protein.interactions import add_distance_threshold
from visualisation.plot import motif_asteroid_plot
from json_conversion.graphein_to_json import g_to_json, load_prot_graph

# Get data 
# --------


#### HARDCODED PROTEIN FOR NOW; TODO: USE INPUT TO ASSIGN PROTEIN ####
#prot_id = "Q9Y2X7"

#prot_id = "4hhb" # PDB file.  Try using with extract_surface_subgraph. 
#g = get_protein_graph(prot_id, config="asa", database)

#psites = get_phosphosites(g)

'''
For protein input sidebar
'''
fnameDict = {'PDB': ['opt1_p', 'opt2_p', 'opt3_p'], 'SWISS_PROT': ['opt1_s', 'opt2_s'], 'AlphaFold': ['opt1_a']}

names = list(fnameDict.keys())
nestedOptions = fnameDict[names[0]]

'''
Generate variables and use get_protein_graph
'''
@callback(
    Output('intermediate-value-prot', 'children'),
    Input('db-dropdown', 'value'),
    Input('prot-input', 'value'),
)
def process_prot_input(db_name, prot_id):
    
    g = get_protein_graph(prot_id, config="asa", database=db_name)
    inter_val = g_to_json(g, prot_id, db_name)
    #return graph in json format
    #print(f"inter_val type is {type(inter_val)}")
    #print(f"inter_val dump type is {type(json.dumps(inter_val))}")
    return json.dumps(inter_val)

'''
Generate psites
'''
@callback(
    Output('intermediate-value-psites', 'data'),
    Input('intermediate-value-prot', 'children'),
    Input('selected-psite-residue-types', 'value'),
    Input('asa-threshold-slider', 'value')
)
def process_psite(graph, sites, rsa_thresh):

    g = load_prot_graph(graph)
    psites = get_phosphosites(g, sites, rsa_thresh)

    return json.dumps(psites)


'''
Radius slider component
'''
def get_marks():

    DEFAULT_MARKING = 10    # 10Å
    BLUE_MARKING    = 3     #  3Å
    
    # Standard markings
    keys = [1, 20, 30]
    vals = [dict(label=f"{k}Å") for k in keys]
    marks = dict(zip(keys, vals))
    
    # Coloured
    marks[DEFAULT_MARKING]  = {'label': f'{DEFAULT_MARKING}Å',  'style': {'color': '#f50'}}
    marks[BLUE_MARKING]     = {'label': f'{BLUE_MARKING}Å',     'style': {'color': '#77b0b1'}}
    
    return marks

'''
ASA Slider component
'''
DEFAULT_ASA_THRESHOLD = 0.5

keys = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
vals = [dict(label=f"{k} ASA") for k in keys]
ASA_THRESHOLD_SLIDER_MARKS = dict(zip(keys, vals))



'''
@callback(Output('intermediate-value', 'data'), Input('selected-psite-dropdown', 'value'))
def update_selected_psite(value):
    
    
    s_g = get_protein_subgraph_radius(g1, site=psite, r=value)
'''


'''
Phosphosite dropdown menu
'''
# @callback(
#     Output('selected-psite-dropdown', 'options'),
#     Input('selected-psite-residue-types', 'value'),
#     Input('intermediate-value-prot', 'children'),
#     )
# def update_psite_dropdown(residues, graph):
# #def update_psite_dropdown(residues):
    
#     tmp = json.loads(graph)
#     g = nx.json_graph.node_link_graph(tmp)
#     g1 = g.copy()
#     return get_phosphosites(g1, residues)
@callback(
    Output('selected-psite-dropdown', 'options'),
    Input('selected-psite-residue-types', 'value'),
    Input('intermediate-value-psites', 'data'),
    )
def update_psite_dropdown(residues, psites):
    tmp = psites.replace('\"', '')
    tmp = tmp.replace('[', '')
    tmp = tmp.replace(']', '')
    options = list(tmp.split(","))
    #print(f"phos dropdown options type is {type(options)}")
    return options

'''
Adjacency matrix plot
'''
@callback(
    Output('graph-adjacency-matrix', 'figure'),
    Output('asteroid-plot', 'figure'),
    Input('radius-threshold-slider', 'value'),
    Input('asa-threshold-slider', 'value'),
    Input('selected-psite-dropdown', 'value'),
    Input('axis-order-dropdown', 'value'),
    Input ('intermediate-value-prot', 'children'),
    )   
    
def update_graphs(radius, asa_threshold, psite, axis_order, graph):
#def update_graph(radius, asa_threshold, psite, axis_order):
    # Get new subgraph
    #print(f"inter psite type is {type(psite)}")

    g1 = load_prot_graph(graph)
    
    if not psite:
        psite = get_phosphosites(g1, residues=['SER','THR','TYR'])[0]
    
    s_g = get_surface_motif(g1, site=psite, r=radius, asa_threshold=asa_threshold)
    # update figure 
    name = g1.graph["name"]
    title = name.upper() + f""" STRUCTURAL MOTIF @ {psite}, threshold: {radius} Å
                                <br>Surface accessibility threshold: {asa_threshold}"""
    adj_mat = get_adjacency_matrix_plot(s_g, psite=psite, title=title, order=axis_order)
    
    # Update asteroid plot
    # Add edges to graph by distance
    #ast_g = add_distance_threshold(s_g, long_interaction_threshold=1, threshold=4.1)
    # Distance edges between AAs closer than 3 A -> so k = radius/3
    k = math.floor(radius/4)
    print (f'Radius = {radius}, {k}-hops shown')
    ast_plt = motif_asteroid_plot(
        g=g1,
        node_id=psite,
        size_nodes_by="rsa",
        node_size_multiplier=80,
        colour_nodes_by="hydrophobicity",
        width=435,
        height=400,
        k=k
    )

    return [adj_mat, ast_plt]


'''
Process Protein Input
'''


'''
Layout
'''
def motifVisualisationTab ():
    return html.Div(className='single-motif-tab-content', children=[
        dcc.Graph(id='graph-adjacency-matrix', className='matrix-display tab-component'),
        dcc.Graph(id='asteroid-plot', className='graph-display tab-component'),
        html.Div(className='options tab-component', children=[
            html.H3('Phosphosite of interest:'),
            dcc.Checklist(id='selected-psite-residue-types',
                options=['SER', 'THR', 'TYR', 'HIS'],
                value=['SER', 'THR', 'TYR'],
                inline=True
            ),
            dcc.Dropdown(
                id='selected-psite-dropdown',
                # options=[{'label':site, 'value':site} for site in psites],
                # value = psites[0],
                style={'width': '80%'}
            ),
            dcc.Dropdown(
                id='axis-order-dropdown',
                options=[{'label':"Sequence position", 'value': "seq"}, 
                    {'label':"Hydrophobicity", 'value':"hydro"}],
                value="hydro",
                style={'width': '80%'}
            )
        ]),
        html.Div(
            id='slider-output-container',
            className='sliders tab-component',
            children=[
                dcc.Slider(0, 30,
                        value=10,
                        marks=get_marks(),
                        included=True, # show trail
                        id='radius-threshold-slider'
                ),
                dcc.Slider(0.0, 1.0,
                        value=DEFAULT_ASA_THRESHOLD,
                        marks=ASA_THRESHOLD_SLIDER_MARKS,
                        included=True, # show trail
                        id='asa-threshold-slider'
                )
            ])
    ])

'''
Debugging to show input
'''
@callback(
    Output('input-show', 'children'),
    Input('db-dropdown', 'value'),
    Input('prot-input', 'value'),
)
def showInput(db, prot_id):
    return u'DB: {} ID: {}'.format(db, prot_id)

def sidebarTab():
    return html.Div(        
        id="sidebar",
        children=[
            html.H3("Find a protein"),
            html.Hr(),
            dcc.Dropdown(
                id='db-dropdown',
                options=[{'label':name, 'value':name} for name in names],
                value = list(fnameDict.keys())[0],
                style={'width': '80%'},
            ),
            # dcc.Dropdown(
            #     id='opt-dropdown',
            #     style={'width': '20%'},
            # ),
            dcc.Input(
                id="prot-input",
                type="text",
                #value="4hhb", # TODO: REMOVE THIS AFTER TESTING
                placeholder="Protein ID",
                style={'width': '80%'},
                debounce=True,
                persistence=True,
                persistence_type='session'
            ),
            dcc.Store(id='intermediate-value-prot', storage_type='session'),
            dcc.Store(id='intermediate-value-psites', storage_type='session'),
            html.Div(id='input-show'), # DEBUGGING inputs
        ],
    )