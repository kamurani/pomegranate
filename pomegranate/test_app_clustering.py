import os
from dash import Dash, html, dcc, Input, Output
import pandas as pd
import plotly.express as px


from definitions import EMBEDDINGS_PATH, STRUCTURE_HUMAN_PATH
from visualisation.plot import motif_plot_distance_matrix
from protein.phosphosite import get_protein_graph

from graphein.protein.visualisation import plot_distance_matrix

import networkx as nx


from typing import Dict, List, Union

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = Dash(__name__, external_stylesheets=external_stylesheets)


"""
Read in data
"""
df = pd.read_csv(EMBEDDINGS_PATH)

graphs: Dict = {str, nx.Graph} # graphs[protein_id]


# TODO: read in graphs from .json savefile 

app.layout = html.Div([
    html.Div([

        html.Div([
            dcc.Dropdown(
                df['Set'].unique(), 
                "Human (known kinases)",
                id='clustering-which-proteome',
            ),
            dcc.RadioItems(
                ['tSNE', 'UMAP'],
                'tSNE',
                id='dim-reduction-method',
                labelStyle={'display': 'inline-block', 'marginTop': '5px'}
            )
        ],
        style={'width': '49%', 'display': 'inline-block'}),

        html.Div([
            dcc.Dropdown(
                df['Method'].unique(),
                'tSNE',
                id='which-visualisation-method'
            ),
            dcc.RadioItems(
                ['TODO', 'TODO'],
                'TODO',
                id='visualisation-options',
                labelStyle={'display': 'inline-block', 'marginTop': '5px'}
            )
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


@app.callback(
    Output('clustering-scatter', 'figure'),
    Input('clustering-which-proteome', 'value'),
    Input('dim-reduction-method', 'value'),
    Input('which-visualisation-method', 'value'),
    Input('visualisation-options', 'value'),
    Input('slider', 'value'))
def update_graph(proteome, dim_reduction_method,
                 visualisation_method, visualisation_option,
                 slider_value):


    # Use slider value to get subset of data

    # TODO:
    # load embeddings appropriate to proteome 
    # run tSNE 

    # FOR NOW: 
    # use pre-calculated clusterings just for hover over. no real time computation yet.  static data for now. 

    xaxis_type = "Linear"
    yaxis_type = "Linear"

    # FILTER
    dff = df[df['Set'] == proteome] 
    dff = dff[dff['Method'] == dim_reduction_method]

    kinase_labels = dff["Kinase"].unique()

    print(kinase_labels)

    fig = px.scatter(
        data_frame=dff,
        x=dff['X'],
        y=dff['Y'],
        hover_name=dff['Protein ID'], #TODO: combine this with psite location to get name on hover.  with name of kinase. 
        custom_data=dff[['Protein ID', 'Phosphosite', 'Kinase']],

    )

    #fig.update_traces(customdata=dff['Kinase'])
    fig.update_traces(customdata=dff['Phosphosite'])

    fig.update_traces(
        hovertemplate="<br>".join([
        "ColX: %{x}",
        "ColY: %{y}",
        "Col1: %{customdata}",
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


@app.callback(
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





if __name__ == '__main__':
    app.run_server(debug=True)