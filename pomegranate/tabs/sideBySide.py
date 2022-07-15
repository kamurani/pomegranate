from dash import Dash, dcc, html, callback, Input, Output

def compareBySideTab (app):
    return html.Div(className='side-by-side-tab-content', children=[
        html.Div(className='first-block tab-component', children=[
            html.Img(src=app.get_asset_url('imgs/gray_thumb.png'), className='first-view'),
            html.Div(className='first-ops', children=[
                html.H3('Options:'),
                dcc.Dropdown(
                    options=[{'label':"Sequence position", 'value': "seq"}, 
                        {'label':"Hydrophobicity", 'value':"hydro"}],
                    value="hydro",
                    style={'width': '80%'}
                )
            ])
        ]),
        html.Div(className='second-block tab-component', children=[
            html.Img(src=app.get_asset_url('imgs/gray_thumb.png'), className='second-view'),
            html.Div(className='second-ops', children=[
                html.H3('Options:'),
                dcc.Dropdown(
                    options=[{'label':"Sequence position", 'value': "seq"}, 
                        {'label':"Hydrophobicity", 'value':"hydro"}],
                    value="hydro",
                    style={'width': '80%'}
                )
            ])
        ])
    ])