"""Web application for pomegranate."""

# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

from distutils.log import debug
import dash
from dash import Dash, dcc, html, Input, Output, State
from tabs.seeMotif import motifVisualisationTab
from tabs.clustering import clustering_tab, update_graph, update_vis_1
from help_tab import help_text
import dash_loading_spinners as dls

import time
#from tabs.sideBySide import compareBySideTab
from tabs.multipleCompare import compareBySideTab

PROTEIN_ID = "default"

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = Dash(__name__, external_stylesheets=external_stylesheets)

'''
Sidebar
'''
fnameDict = {'PDB': ['opt1_p', 'opt2_p', 'opt3_p'], 'SWISS_PROT': ['opt1_s', 'opt2_s'], 'AlphaFold': ['opt1_a']}

names = list(fnameDict.keys())
nestedOptions = fnameDict[names[0]]

sidebar = html.Div(
    id="sidebar-container",
    children=[
        html.H3("Find a protein"),
        html.Hr(),
        dcc.Dropdown(
            id='db-dropdown',
            options=[{'label':name, 'value':name} for name in names],
            value = list(fnameDict.keys())[0],
            style={'width': '80%'},
        ),
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

'''
Tabs
'''
tab_selected_style = {
    'borderTop': '3px solid #b52d37',
}

content = html.Div(
    children= [
        dcc.Tabs(id="tab-options", value='single-motif-view', children=[
            dcc.Tab(label='Visualisation', value='single-motif-view', selected_style=tab_selected_style),
            dcc.Tab(label='Compare Motifs', value='multi-motif-view', selected_style=tab_selected_style),
            dcc.Tab(label='Proteome View', value='clustering', selected_style=tab_selected_style),
            dcc.Tab(label='Help', value='documentation', selected_style=tab_selected_style),
        ]),
        html.Div(id='tab-container'),
    ]
)

base_page = html.Div(
    children=[
        html.Div(
            id="div-loading",
            children=[
                dls.Pacman(
                    fullscreen=True, 
                    id="loading-whole-app"
                )
            ]
        ),
        html.Div(
            className="div-app",
            id="div-app",
            children = [ html.Div([
                    html.Img(src=app.get_asset_url('imgs/POMEGRANATE-LOGO.png'), style={'width': '40%'}),
                    html.H2('PhOsphosite Motif Explorer -- GRAph Network Abstraction Through Embeddings'),
                    html.Div(id="content-grid", children=[sidebar,content]),
                ])
        ]
        )
    ]
)

app.layout = base_page

@app.callback(
    Output("div-loading", "children"),
    [
        Input("div-app", "loading_state")
    ],
    [
        State("div-loading", "children"),
    ]
)
def hide_loading_after_startup(
    loading_state, 
    children
    ):
    if children:
        print("remove loading spinner!")
        return None
    print("spinner already gone!")
    raise PreventUpdate


app.validation_layout = html.Div([
    base_page,
    motifVisualisationTab()
])

@app.callback(Output('tab-container', 'children'),
              Input('tab-options', 'value'))
def render_content(tab):
    if tab == 'single-motif-view':
        time.sleep(2);
        return motifVisualisationTab()
    elif tab == 'multi-motif-view':
        return compareBySideTab()
    elif tab == 'clustering':
        return clustering_tab()
    elif tab == 'documentation':
        return help_text()




def run():
    app.run_server(debug=True, port=8050)

if __name__ == '__main__':
    run()
