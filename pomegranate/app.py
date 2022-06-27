"""Web application for pomegranate."""

# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import dash
from dash import Dash, dcc, html

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = Dash(__name__, external_stylesheets=external_stylesheets, use_pages=True)

app.layout = html.Div([
	html.H1('POMEGRANATE'),
    html.H2('PhOsphosite Motif Explorer -- GRAph Network Abstraction Through Embeddings')

    html.Div(
        [
            html.Div(
                dcc.Link(
                    f"{page['name']} - {page['path']}", href=page["relative_path"]
                )
            )
            for page in dash.page_registry.values()
        ]
    ),

	dash.page_container
])



if __name__ == '__main__':
    app.run_server(debug=True)
