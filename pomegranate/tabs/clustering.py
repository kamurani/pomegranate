"""Clustering tab"""
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

import os
from xml.etree.ElementInclude import include
from dash import Dash, html, dcc, Input, Output
import pandas as pd
import plotly.express as px


from definitions import EMBEDDINGS_PATH, STRUCTURE_HUMAN_PATH
from utils.amino_acid import aa1letter
from visualisation.plot import motif_plot_distance_matrix
from protein.phosphosite import get_protein_graph

from graphein.protein.visualisation import plot_distance_matrix
from graphein.utils.utils import protein_letters_3to1_all_caps as aa3to1


import networkx as nx


from typing import Callable, Dict, List, Union


"""
Read in data
"""
df = pd.read_csv(EMBEDDINGS_PATH)

graphs: Dict[str, nx.Graph] # graphs[protein_id]


"""
Layout
"""
def clustering_tab():
    return html.Div([
        html.Div([
            
            # row
            html.Div([
                html.Div(
                    [
                        html.Label(
                            ['Method'], style={'font-weight': 'bold', "text-align": "left"}
                        ),
                        dcc.RadioItems(
                            ['tSNE', 'UMAP'],
                            'tSNE',
                            id='dim-reduction-method',
                            labelStyle={'display': 'inline-block', 'marginTop': '5px'}
                        ),
                    ], style={'display': 'inline-block'},
                ),
                html.Div(
                    [
                        html.Label(
                            ['Phosphosite Residues'], style={'font-weight': 'bold', "text-align": "left"}
                        ),
                        dcc.Checklist(id='selected-psite-residue-types',
                            options=['SER', 'THR', 'TYR', 'HIS'],
                            value=['SER', 'THR', 'TYR'],
                            inline=True,
                            style={'display': 'inline-block'},
                            labelStyle={'display': 'inline-block', 'marginTop': '5px'},
                        ),
                    ], style={"margin-left": "15px", 'display': 'inline-block'},
                ),
                html.Div(
                    [
                        html.Label(
                            ['Colour by'], style={'font-weight': 'bold', "text-align": "left"}
                        ),                   
                        
                        dcc.RadioItems(id='colour-by',
                            options=['Residue', 'Kinase', 'Average RSA', 'Cluster'],
                            value='Kinase',
                            inline=True,
                            labelStyle={'display': 'inline-block', 'marginTop': '5px'},
                            style={'display': 'inline-block'}
                        ),
                    ], style={"margin-left": "20px", 'display': 'inline-block'},
                ),


            ]),
            # row 
            # TODO

            

            html.Div([
                html.Label(
                    ['Proteome'], style={'font-weight': 'bold', "text-align": "left"}
                ),
                dcc.Dropdown(
                    df['Set'].unique(), 
                    "Human (known kinases)",
                    id='clustering-which-proteome',
                ),
                
            ],
            style={'width': '49%', 'display': 'inline-block'}),

            html.Div([
                dcc.Dropdown(
                    df['Method'].unique(),
                    'tSNE',
                    id='which-visualisation-method'
                ),
                
            ], style={'width': '49%', 'float': 'right', 'display': 'inline-block'})
        ], style={
            'padding': '10px 5px'
        }),

        # Clustering plot
        html.Div([
            dcc.Graph(
                id='clustering-scatter',
                hoverData={'points': [{'customdata': 'DEFAULT', 'psite': 'DEFAULT'}]} # TODO
            )
        ], style={'width': '49%', 'display': 'inline-block', 'padding': '0 20'}),

        
        # Visualisation plots
        html.Div([
            dcc.Graph(id='vis-1'),
            #dcc.Graph(id='vis-2'), # TODO, 2nd visualisation graph
        ], style={'display': 'inline-block', 'width': '49%'}),

        

        # Slider TODO
        html.Div(dcc.Slider(
            0,
            10,
            step=None,
            id='slider',
            value=10,
        ), style={'width': '49%', 'padding': '0px 20px 20px 20px'})
    ])


@callback(
    Output('clustering-scatter', 'figure'),
    Input('selected-psite-residue-types', 'value'),
    Input('clustering-which-proteome', 'value'),
    Input('dim-reduction-method', 'value'),
    Input('which-visualisation-method', 'value'),
    Input('colour-by', 'value'),
    Input('slider', 'value'),
)
def update_graph(include_residues, proteome, dim_reduction_method,
                 visualisation_method, colour_by,
                 slider_value):

    print("Residues:",include_residues)

    residues = include_residues
    # Use slider value to get subset of data

    # TODO:
    # load embeddings appropriate to proteome 
    # run tSNE 

    # FOR NOW: 
    # use pre-calculated clusterings just for hover over. no real time computation yet.  static data for now. 

    xaxis_type = "Linear"
    yaxis_type = "Linear"

    # SPLIT DATAFRAME 
    df[['Chain ID', 
        'Phosphosite Residue',
        'Phosphosite Sequence Position']] = df.Phosphosite.apply(
            lambda x: pd.Series(str(x).split(':'))
        )

    # FILTER
    dff = df[df['Set'] == proteome] 


    

    """Returns a function that can filter a dataframe of residue IDs"""
    def get_residue_filter(
        residues: Union[List[str], str],
        invert: bool = False,
    ) -> Callable:
        """
        :param residues: Either a string containing 1-letter codes. 
        :type residues: Union[List[str], str] 
        :param invert: Return True if the input is NOT in the specified list of residues.  Defaults to ``False``. 
        :type invert: bool


        To allow for all residues (i.e. apply no filtering), an empty list can be supplied with ``invert`` set to ``True``.
        
        """
        if type(residues) == str: # String containing 1-letter codes
            residues = residues.upper()
        else:   # List
            residues = "".join(aa1letter(x.upper()) for x in residues)


        return lambda x: not invert if (aa3to1(x.split(':')[1]) in residues) else invert

        A, B = True, False 
        if invert: A, B = B, A


    # If no residues are selected, include all. 
    filt = get_residue_filter(residues, invert=(not residues))

    
    if colour_by == "Residue": colour = 'Phosphosite Residue'
    else: colour = colour_by

    dff = dff[dff['Phosphosite'].apply(filt)]
    dff = dff[dff['Method'] == dim_reduction_method]

    kinase_labels = dff["Kinase"].unique()

    print(kinase_labels)

    fig = px.scatter(
        #name=f"{dim_reduction_method} clustering of motifs from {proteome}",
        data_frame=dff,
        x='X',
        y='Y',
        color=colour,

        hover_name=dff['Protein ID'], #TODO: combine this with psite location to get name on hover.  with name of kinase. 


        #custom_data=dff[['Protein ID', 'Phosphosite', 'Kinase']],

    )

    #fig.update_traces(customdata=dff['Kinase'])
    fig.update_traces(customdata=dff[['Protein ID', 'Phosphosite', 'Kinase']])

    fig.update_traces(
        hovertemplate="<br>".join([
        "%{customdata[0]}",
        "Kinase: %{customdata[2]}",
        "",
        "X: %{x}",
        "Y: %{y}",
        
        #"Col2: %{customdata[1]}",
    ])
    )

    x_axis, y_axis = dim_reduction_method+"1", dim_reduction_method+"2"
    fig.update_xaxes(title=x_axis, type='linear' if xaxis_type == 'Linear' else 'log')
    fig.update_yaxes(title=y_axis, type='linear' if yaxis_type == 'Linear' else 'log')

    fig.update_layout(margin={'l': 40, 'b': 40, 't': 10, 'r': 0}, hovermode='closest')
    

    return fig


def create_time_series(dff, axis_type, title):

    fig = px.scatter(dff, x='Year', y='Value')

    fig.update_traces(mode='lines+markers')

    fig.update_xaxes(showgrid=False)

    fig.update_yaxes(type='linear' if axis_type == 'Linear' else 'log')

    fig.add_annotation(x=0, y=0.85, xanchor='left', yanchor='bottom',
                       xref='paper', yref='paper', showarrow=False, align='left',
                       text=title)

    fig.update_layout(height=225, margin={'l': 20, 'b': 30, 'r': 10, 't': 10})

    return fig


@callback(
    Output('vis-1', 'figure'),
    Input('clustering-scatter', 'hoverData'),
    Input('which-visualisation-method', 'value'),
    Input('visualisation-options', 'value'),

)
def update_vis_1(
    hoverData,
    vis_method,
    vis_options, 

):

    # Load graph
    # may have to do this from file as we go if not enough memory? 
    #g = graphs[]

    # TODO: present kinase name and psite location on hover_over



    data = hoverData['points'][0]['customdata']

    print(f"data: {data}")


    name = data[0]
    if name is None or name == "DEFAULT": 
        fig = {}
        return fig

    protein_id = name.split('@')[0].strip()


    try: 
        g = graphs[protein_id]
    except:
        pdb_path = os.path.join(STRUCTURE_HUMAN_PATH, f"{protein_id}.pdb")
        print(f"Constructing graph from {pdb_path}...")
        g = get_protein_graph(config="rsa", pdb_path=pdb_path)
    
    try: 
        psite = data[1]
    except:
        psite = None

    #if psite: fig = motif_plot_distance_matrix(g=g, psite=psite, aa_order="hydro", show_residue_labels=False)
    #else: fig = plot_distance_matrix(g)
    # Subgraph 

    #psite = 


    # Get figure
    #fig = motif_plot_distance_matrix(g=g, psite=psite)

    
    
    return fig