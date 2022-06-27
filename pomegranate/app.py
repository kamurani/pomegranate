"""Web application for pomegranate."""

# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

from dash import Dash, dcc, html, Input, Output
import plotly.express as px
import pandas as pd
import numpy as np

from protein.phosphosite import *   





external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = Dash(__name__, external_stylesheets=external_stylesheets)


# Get data 
# --------

prot_id = "Q9Y2X7"
g = get_protein_graph(prot_id)

psites = get_phosphosites(g)



'''
Radius slider component
'''
def get_marks():

    BLUE_MARKING    = 3     #  3Å
    DEFAULT_MARKING = 10    # 10Å
    
    # Standard markings
    keys = [1, 20, 30]
    vals = [dict(label=f"{k}Å") for k in keys]
    marks = dict(zip(keys, vals))
    
    # Coloured
    marks[DEFAULT_MARKING]  = {'label': f'{DEFAULT_MARKING}Å',  'style': {'color': '#f50'}}
    marks[BLUE_MARKING]     = {'label': f'{BLUE_MARKING}Å',     'style': {'color': '#77b0b1'}}
    
    return marks

'''
Layout
'''
app.layout = html.Div([
    dcc.Graph(id='graph-adjacency-matrix'),
    dcc.Slider(0, 30,
               value=10,
               marks=get_marks(),
               included=True, # show trail
               id='radius-threshold-slider'
    ),
    
    html.Div([
        dcc.Checklist(id='selected-psite-residue-types',
                    options=['SER', 'THR', 'TYR', 'HIS'],
                    value=['SER', 'THR', 'TYR'],
                    inline=True
        ),
        dcc.Dropdown(
            id='selected-psite-dropdown',
            options=[{'label':site, 'value':site} for site in psites],
            value = psites[0]
            ),
    ], style={'width': '20%', 'display': 'inline-block'}
    ),
            
    html.Div([
    dcc.Dropdown(
        id='opt-dropdown',
        ),
        ],style={'width': '20%', 'display': 'inline-block'}
    ),
    
    html.Div(id='slider-output-container')
])





'''
@app.callback(Output('intermediate-value', 'data'), Input('selected-psite-dropdown', 'value'))
def update_selected_psite(value):
    
    
    s_g = get_protein_subgraph_radius(g1, site=psite, r=value)
'''

@app.callback(
    Output('selected-psite-dropdown', 'options'),
    Input('selected-psite-residue-types', 'value'),
    )
def update_psite_dropdown(residues):
    
    g1 = g.copy()
    return get_phosphosites(g1, residues)


@app.callback(
    Output('graph-adjacency-matrix', 'figure'),
    Input('radius-threshold-slider', 'value'),
    Input('selected-psite-dropdown', 'value'),
    )
    
def update_graph(radius, psite):
    # Get new subgraph
    g1 = g.copy()
    
    if not psite:
        psite = get_phosphosites(g1)[0]
    
    s_g = get_protein_subgraph_radius(g1, site=psite, r=radius)
    # update figure 
    title = g.graph["name"] + f" MOTIF @ {psite}, threshold = {radius} Å"
    figure = get_adjacency_matrix_plot(s_g, psite=psite, title=title)
    
    
    return figure
    



if __name__ == '__main__':
    app.run_server(debug=True)
