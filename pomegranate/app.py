"""Web application for pomegranate."""

# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import dash
from dash import Dash, dcc, html, Input, Output
from tabs.seeMotif import motifVisualisationTab, sidebarTab

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = Dash(__name__, external_stylesheets=external_stylesheets)

# sidebar = html.Div(
#     id="sidebar",
#     children= [
#         html.H3("Find a protein"),
#         html.Hr(),
#         html.P(
#             "A placeholder for searching a protein of interest"
#         ),
#         html.Button("Button")
#     ]
# )

sidebar = html.Div(
    id="sidebar-container",
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
        return html.H3('Documentation')

'''
NOT SURE WHAT TO DO FOR INPUT
'''
@app.callback(Output('sidebar-container', 'children'),
              Input('tab-options', 'value'))
def render_content(tab):
    return sidebarTab()


if __name__ == '__main__':
    app.run_server(debug=True)
