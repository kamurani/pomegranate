from dash import Dash, dcc, html, callback, Input, Output

from protein.phosphosite import get_adjacency_matrix_plot, get_phosphosites, get_protein_graph, get_surface_motif

from visualisation.plot import multiple_motif_plot_distance_matrix

# Get data 
# --------
#### HARDCODED PROTEIN FOR NOW; TODO: USE INPUT TO ASSIGN PROTEIN ####
prot_id = "4hhb" # PDB file.  Try using with extract_surface_subgraph. 

USE_ALPHAFOLD = False
g = get_protein_graph(prot_id, use_alphafold=USE_ALPHAFOLD, config="asa")
psites = get_phosphosites(g)


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

'''
ASA Slider component
'''
DEFAULT_ASA_THRESHOLD = 0.5

keys = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
vals = [dict(label=f"{k} ASA") for k in keys]
ASA_THRESHOLD_SLIDER_MARKS = dict(zip(keys, vals))
#####################################################################

'''
Adjacency matrix plot
'''
@callback(
    Output('adj-matrices', 'figure'),
    Input('radius-threshold-slider', 'value'),
    Input('asa-threshold-slider', 'value'),
    Input('axis-order-dropdown', 'value'),
    Input('colour-dropdown', 'value')
    )
    
def update_graph(radius, asa_threshold, axis_order, colour):
    
    # Get new subgraph
    g1 = g.copy()
    # Get phosphosite
    # TODO get phosphosite info. Probably want to select this before comparing?
    psite = get_phosphosites(g1)[0]
    
    s_g = get_surface_motif(g1, site=psite, r=radius, asa_threshold=asa_threshold)
    # update figure 
    name = g.graph["name"]
    title = name.upper() + f""" STRUCTURAL MOTIF @ {psite}"""
    figure = multiple_motif_plot_distance_matrix(to_plot=[(s_g, psite), (s_g, psite), (s_g, psite), (s_g, psite)], colour=colour)
      
    return figure


'''
Layout
'''
def compareBySideTab ():
    return html.Div(className='side-by-side-tab-content', children=[
        html.Div(className='tab-component', children=[
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
            ),
            html.Div(children=[
                html.H4('Order matrix by ', style={'float':'left', 'margin':'5px', 'height':'20px', 'bottom':'0'}),
                dcc.Dropdown(
                    id='axis-order-dropdown',
                    options=[{'label':"Sequence position", 'value': "seq"}, 
                        {'label':"Hydrophobicity", 'value':"hydro"}],
                    value="hydro",
                    style={'float':'left', 'margin':'5px', 'height':'20px', 'bottom':'0', 'width':'200px'}
                ),
                html.H4('Grayscale? ', style={'float':'left', 'margin':'5px', 'height':'20px', 'bottom':'0'}),
                dcc.Dropdown(
                    id='colour-dropdown',
                    options=[{'label':"No", 'value': "viridis_r"}, 
                        {'label':"Yes", 'value':"gray"}],
                    value="viridis_r",
                    style={'float':'left', 'margin':'5px', 'height':'20px', 'bottom':'0', 'width':'200px'}
                )],
                style={'height':'30px'}
            )
        ]),
        dcc.Graph(id='adj-matrices', className='tab-component')
    ])