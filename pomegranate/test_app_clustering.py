from dash import Dash, html, dcc, Input, Output
import pandas as pd
import plotly.express as px


from definitions import EMBEDDINGS_PATH

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = Dash(__name__, external_stylesheets=external_stylesheets)


"""
Read in data
"""
df = pd.read_csv(EMBEDDINGS_PATH)


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
            hoverData={'points': [{'customdata': 'Japan'}]} # TODO
        )
    ], style={'width': '49%', 'display': 'inline-block', 'padding': '0 20'}),

    # Visualisation plots
    html.Div([
        dcc.Graph(id='vis-1'),
        dcc.Graph(id='vis-2'),
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
    dff = df[df['Method'] == dim_reduction_method]

    kinase_labels = dff[dff["Kinase"]].unique()

    print(kinase_labels)
    exit

    fig = px.scatter(
        x=dff['X'],
        y=dff['Y'],
        hover_name=dff['Protein ID'], #TODO: combine this with psite location to get name on hover.  with name of kinase. 
    )

    """
    fig.update_traces(customdata=dff[dff['Indicator Name'] == yaxis_column_name]['Country Name'])

    fig.update_xaxes(title=xaxis_column_name, type='linear' if xaxis_type == 'Linear' else 'log')

    fig.update_yaxes(title=yaxis_column_name, type='linear' if yaxis_type == 'Linear' else 'log')

    fig.update_layout(margin={'l': 40, 'b': 40, 't': 10, 'r': 0}, hovermode='closest')
    """

    return fig


"""

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
    Output('x-time-series', 'figure'),
    Input('crossfilter-indicator-scatter', 'hoverData'),
    Input('crossfilter-xaxis-column', 'value'),
    Input('crossfilter-xaxis-type', 'value'))
def update_y_timeseries(hoverData, xaxis_column_name, axis_type):
    country_name = hoverData['points'][0]['customdata']
    dff = df[df['Country Name'] == country_name]
    dff = dff[dff['Indicator Name'] == xaxis_column_name]
    title = '<b>{}</b><br>{}'.format(country_name, xaxis_column_name)
    return create_time_series(dff, axis_type, title)


@app.callback(
    Output('y-time-series', 'figure'),
    Input('crossfilter-indicator-scatter', 'hoverData'),
    Input('crossfilter-yaxis-column', 'value'),
    Input('crossfilter-yaxis-type', 'value'))
def update_x_timeseries(hoverData, yaxis_column_name, axis_type):
    dff = df[df['Country Name'] == hoverData['points'][0]['customdata']]
    dff = dff[dff['Indicator Name'] == yaxis_column_name]
    return create_time_series(dff, axis_type, yaxis_column_name)

"""

if __name__ == '__main__':
    app.run_server(debug=True)