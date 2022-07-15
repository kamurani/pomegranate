"""Web application for pomegranate."""

# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

from distutils.log import debug
import dash
from dash import Dash, dcc, html, Input, Output
from tabs.seeMotif import motifVisualisationTab
from tabs.sideBySide import compareBySideTab

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

app.layout = html.Div([
    html.Img(src=app.get_asset_url('imgs/POMEGRANATE-LOGO.png'), style={'width': '40%'}),
    html.H2('PhOsphosite Motif Explorer -- GRAph Network Abstraction Through Embeddings'),
    html.Div(id="content-grid", children=[sidebar,content])
    
])

@app.callback(Output('tab-container', 'children'),
              Input('tab-options', 'value'))
def render_content(tab):
    if tab == 'single-motif-view':
        return motifVisualisationTab()
    elif tab == 'multi-motif-view':
        return compareBySideTab()
    elif tab == 'clustering':
        return html.H3('Look how cool our clusters are')
    elif tab == 'documentation':
        return html.H3('Documentation')



def run():
    app.run_server(debug=True, port=8050)

if __name__ == '__main__':
    run()
