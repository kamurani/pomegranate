from pydoc import classname
from dash import Dash, dcc, html, callback, Input, Output
import math

from json_conversion.graphein_to_json import load_prot_graph
from protein.phosphosite import get_phosphosites, get_surface_motif
from visualisation.plot import multiple_motif_plot_distance_matrix, multi_asteroid_plot


# Copied from seeMotif ... is there a better way to hold this info? #
#####################################################################
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
    marks[DEFAULT_MARKING]  = {'label': f'{DEFAULT_MARKING}Å',  'style': {'color': '#b52d37'}}
    marks[BLUE_MARKING]     = {'label': f'{BLUE_MARKING}Å',     'style': {'color': '#77b0b1'}}
    
    return marks
#######################################################################

@callback(
    Output('psites-to-compare', 'options'),
    Output('psites-to-compare', 'value'),
    Input('intermediate-value-prot', 'children'),
    )
def update_offered_psites(graph):

    g = load_prot_graph(graph)
    psites = get_phosphosites(g)

    return psites, [psites[0]]


'''
Adjacency matrix plot
'''
@callback(
    Output('adj-matrices', 'figure'),
    Input('radius-threshold-slider', 'value'),
    Input('axis-order-dropdown-cf', 'value'),
    Input('intermediate-value-prot', 'children'),
    Input('psites-to-compare', 'value'),
    Input('colour-dropdown', 'value'),
    Input('cf-vis-method', 'value')
    )  
def update_graph(radius, axis_order, graph, psites, colour, vis_type):
    
    # Get new subgraph
    g = load_prot_graph(graph)
    
    # Get phosphosite
    # TODO get phosphosite info. Probably want to select this before comparing?
    graphs_to_plot = []
    for psite in psites:
        s_g = get_surface_motif(g.copy(), site=psite, r=radius, asa_threshold=0.5)
        graphs_to_plot.append((s_g, psite))

    # update figure
    if vis_type == 'matrix':
        figure = multiple_motif_plot_distance_matrix(to_plot=graphs_to_plot, aa_order=axis_order, colour=colour)
    elif vis_type == 'asteroid':
        k = math.floor(radius/4)
        figure = multi_asteroid_plot(to_plot=graphs_to_plot, k=k, colour=colour)
      
    return figure


'''
Layout
'''
def compareBySideTab ():
    return html.Div(className='side-by-side-tab-content', children=[
        html.Div(className='options-grid', children=[
            html.Div(className='select-multi-psites tab-component', children=[
                html.H3('Phosphosites to compare: '),
                dcc.Dropdown(
                    id='psites-to-compare',
                    multi=True,
                )
            ]),
            html.Table(className='cf-customising tab-component', children=[
                html.Tr([
                    html.Td([html.H5('Plot as ')], style={'text-align':'right'}),
                    html.Td([dcc.Dropdown(
                        id='cf-vis-method',
                        options=[{'label':'Adjacency Matrix', 'value':'matrix'},
                                {'label':'Asteroid Plot', 'value':'asteroid'}],
                        value='matrix',
                        style={'height':'20px', 'bottom':'0', 'width':'200px'}
                    )]),
                    html.Td([html.H5('Grayscale? ')], style={'text-align':'right'}),
                    html.Td([dcc.Dropdown(
                        id='colour-dropdown',
                        options=[{'label':"No", 'value': "amp"}, 
                            {'label':"Yes", 'value':"gray_r"}],
                        value="amp",
                        style={'height':'20px', 'bottom':'0', 'width':'200px'}
                    )]),
                    html.Td([html.H5('Order matrix by ')], style={'text-align':'right'}),
                    html.Td([dcc.Dropdown(
                        id='axis-order-dropdown-cf',
                        options=[{'label':"Sequence position", 'value': "seq"}, 
                            {'label':"Hydrophobicity", 'value':"hydro"}],
                        value="hydro",
                        style={'height':'20px', 'width':'200px'}
                    )])
                ])
            ]),
            html.Div(className='cf-radius-slider tab-component', children=[
                html.Br(),
                dcc.Slider(0, 30,
                    value=10,
                    marks=get_marks(),
                    included=True, # show trail
                    id='radius-threshold-slider',
            )
            ])
        ]),
        html.Div(className='tab-component', children=[
            dcc.Graph(id='adj-matrices')
        ])       
    ])