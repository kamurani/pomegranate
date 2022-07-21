"""Web application for pomegranate."""

# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

from distutils.log import debug
import dash
from dash import Dash, dcc, html, Input, Output, State
from tabs.seeMotif import motifVisualisationTab
from help_tab import help_text

import dash_loading_spinners as dls

import time

PROTEIN_ID = "default"

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = Dash(__name__, external_stylesheets=external_stylesheets)

sidebar = html.Div(
    id="sidebar",
    children= [
        html.H3("Find a protein"),
        html.Hr(),
        html.P(
            "A placeholder for searching a protein of interest"
        ),
        html.Button("Button")
    ]
)

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

app.layout = html.Div(
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

@app.callback(Output('tab-container', 'children'), Input('tab-options', 'value'))
def render_content(tab):
    if tab == 'single-motif-view':
        time.sleep(2);
        return motifVisualisationTab()
    elif tab == 'multi-motif-view':
        return html.Div([
            html.H3('Placeholder graph'),
            dcc.Graph(
                figure={
                    'data': [{
                        'x': [1, 2, 3],
                        'y': [3, 1, 2],
                        'type': 'bar'
                    }]
                }
            )
        ])
    elif tab == 'clustering':
        return html.H3('Look how cool our clusters are')
    elif tab == 'documentation':
        return help_text()



def run():
    app.run_server(debug=True, port=8050)

if __name__ == '__main__':
    run()
