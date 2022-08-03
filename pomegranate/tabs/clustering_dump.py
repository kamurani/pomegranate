"""Clustering tab"""
from dash import Dash, dcc, html, callback, Input, Output
import plotly.express as px
import pandas as pd
import numpy as np

from protein.phosphosite import *
from protein.phosphosite import get_surface_motif  



'''
@callback(Output('intermediate-value', 'data'), Input('selected-psite-dropdown', 'value'))
def update_selected_psite(value):
    
    
    s_g = get_protein_subgraph_radius(g1, site=psite, r=value)
'''


'''
Phosphosite dropdown menu
'''
@callback(
    Output('selected-psite-dropdown', 'options'),
    Input('selected-psite-residue-types', 'value'),
    )
def update_psite_dropdown(residues):
    
    g1 = g.copy()
    return get_phosphosites(g1, residues)


'''
Adjacency matrix plot
'''
@callback(
    Output('graph-adjacency-matrix', 'figure'),
    Input('radius-threshold-slider', 'value'),
    Input('asa-threshold-slider', 'value'),
    Input('selected-psite-dropdown', 'value'),
    Input('axis-order-dropdown', 'value'),
    )
    
def update_graph(radius, asa_threshold, psite, axis_order):
    # Get new subgraph
    g1 = g.copy()
    
    if not psite:
        psite = get_phosphosites(g1)[0]
    
    s_g = get_surface_motif(g1, site=psite, r=radius, asa_threshold=asa_threshold)
    # update figure 
    name = g.graph["name"]
    title = name.upper() + f""" STRUCTURAL MOTIF @ {psite}, threshold: {radius} Å
                                <br>Surface accessibility threshold: {asa_threshold}"""
    figure = get_adjacency_matrix_plot(s_g, psite=psite, title=title, order=axis_order)
      
    return figure

'''
Layout
'''
def motifVisualisationTab ():
    return html.Div(className='single-motif-tab-content', children=[
        dcc.Graph(id='graph-adjacency-matrix', className='matrix-display tab-component'),
        html.Div(className='graph-display tab-component', children=[
            html.H4('show graph with surface mesh here')
        ]),
        html.Div(className='options tab-component', children=[
            html.H3('Phosphosite of interest:'),
            dcc.Checklist(id='selected-psite-residue-types',
                options=['SER', 'THR', 'TYR', 'HIS'],
                value=['SER', 'THR', 'TYR'],
                inline=True
            ),
            dcc.Dropdown(
                id='selected-psite-dropdown',
                options=[{'label':site, 'value':site} for site in psites],
                value = psites[0],
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