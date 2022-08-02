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
from visualisation.plot import asteroid_plot_2, motif_asteroid_plot
from json_conversion.graphein_to_json import g_to_json, load_prot_graph

import os
from xml.etree.ElementInclude import include
from dash import Dash, html, dcc, Input, Output
import pandas as pd
import plotly.express as px


from definitions import EMBEDDINGS_PATH, STRUCTURE_HUMAN_PATH, SAVED_CLUSTER_GRAPHS_PATH, STRUCTURE_YEAST_PATH
from utils.amino_acid import aa1letter
from visualisation.plot import motif_plot_distance_matrix, motif_asteroid_plot
from protein.phosphosite import get_protein_graph

from graphein.protein.visualisation import plot_distance_matrix
from graphein.utils.utils import protein_letters_3to1_all_caps as aa3to1

from utils.clustering_functions import construct_graphs


import networkx as nx


from typing import Callable, Dict, List, Union


"""
Read in data
"""
df = pd.read_csv(EMBEDDINGS_PATH)



graphs: Dict[str, nx.Graph]


graphs = construct_graphs(
    df=df,
    pdb_dir=STRUCTURE_HUMAN_PATH,
    out_dir=SAVED_CLUSTER_GRAPHS_PATH,
    overwrite=False,
) 

graphs = construct_graphs(
    df=df,
    pdb_dir=STRUCTURE_YEAST_PATH,
    out_dir=SAVED_CLUSTER_GRAPHS_PATH,
    overwrite=False,
) 
# graphs[protein_id]





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
                            options=['Residue', 'Kinase', 'Kinase (discrete)'],
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
                html.Label(
                    ['Visualisation'], style={'font-weight': 'bold', "text-align": "left"}
                ),
                dcc.Dropdown(
                    ['Matrix', 'Asteroid'],
                    'Matrix',
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
        html.Div([
            
            html.Label(
                    ['Radius threshold'], style={'font-weight': 'bold', "text-align": "left"}
                ),

            dcc.Slider(
                0,
                10,
                step=None,
                id='slider',
                value=10,
            ), 
            ],
            style={'width': '40%', 'padding': '0px 20px 20px 20px'})
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

    #print("Residues:",include_residues)

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
    
    dff = dff[dff['Method'] == dim_reduction_method]


    kin_to_num: Dict[str, int] = {}
    for i, k in enumerate(dff['Kinase'].unique()):
        kin_to_num[k] = i


    dff['Kinase Number'] = dff.Kinase.apply(
        lambda x: kin_to_num[x]
    )



    

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
    elif colour_by == "Kinase": colour = "Kinase Number"
    elif colour_by == "Kinase (discrete)": colour = "Kinase"
    else: colour = colour_by

    dff = dff[dff['Phosphosite'].apply(filt)]
    dff = dff[dff['Method'] == dim_reduction_method]

    #kinase_labels = dff["Kinase"].unique()

    #print(kinase_labels)
    #print(f"{len(kinase_labels)} kinases.")

    fig = px.scatter(
        #name=f"{dim_reduction_method} clustering of motifs from {proteome}",
        data_frame=dff,
        x='X',
        y='Y',
        color=colour,
        color_continuous_scale=px.colors.sequential.Viridis,
        color_discrete_sequence=px.colors.sequential.Jet,

        hover_name=dff['Protein ID'], #TODO: combine this with psite location to get name on hover.  with name of kinase. 


        custom_data=dff[['Protein ID', 'Phosphosite', 'Kinase']],

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
    Input('slider', 'value'),
    Input('clustering-which-proteome', 'value'),

)
def update_vis_1(
    hoverData,
    vis_method, 
    radius,
    proteome,
):

    # Load graph
    # may have to do this from file as we go if not enough memory? 
    #g = graphs[]

    # TODO: present kinase name and psite location on hover_over

    data = hoverData['points'][0]['customdata']
    print(f"data: {data}")

    name        = data[0]
    phosphosite = data[1]
    kinase      = data[2]

    if name is None or name == "DEFAULT": 
        fig = {}
        return fig

    protein_id = name.split('@')[0].strip()

    database = "AF2"
    prot_name = f"{protein_id}_{phosphosite}"
    
    filename = f"{prot_name}_{database}.json"
    prot_save_path = os.path.join(SAVED_CLUSTER_GRAPHS_PATH, filename)

    # Check if exists 
    if not os.path.isfile(prot_save_path): return {}

    
    print(f"loading from {prot_save_path} ...")
    with open(prot_save_path, 'r') as f:
        data = json.load(f)
    
    print(type(data))
    g: nx.Graph = load_prot_graph(data)

    


    
    
    if True:
        if vis_method == 'Matrix':
            g = get_protein_subgraph_radius(g=g, site=phosphosite, r=radius)
            fig = motif_plot_distance_matrix(g=g, psite=phosphosite, aa_order="hydro", show_residue_labels=True)
        elif vis_method == 'Asteroid':
            print('Asteroid plot.')
            k = min(1, math.floor(radius/4))
            fig = motif_asteroid_plot(
                g=g,
                
                node_id=phosphosite,
                size_nodes_by="rsa",
                node_size_multiplier=80,
                colour_nodes_by="hydrophobicity",
                #width=435,
                #height=400,
                k=k,
            )
    try:
        pass
    except:
        fig = {}
    #if psite: fig = motif_plot_distance_matrix(g=g, psite=psite, aa_order="hydro", show_residue_labels=False)
    #else: fig = plot_distance_matrix(g)
    # Subgraph 

    #psite = 


    # Get figure
    #fig = motif_plot_distance_matrix(g=g, psite=psite)


    return fig