import dash
import dash_bootstrap_components as dbc
from dash import html, dcc
from dash.dependencies import Input, Output
#connect to main app.py file
from app import app
from app import server

app.title = "CIMEA-VEC Reports"

#connect to your app pages
from apps import home,sequencing, global_updates,run231,run232,run233,run234,run235,run237,run238,run239

runs = ["230", "231", "232", "233", "234", "235", "236", "237", "238",'239']

dropdown_items = [dbc.DropdownMenuItem(f"Run{run}", href=f"/apps/sequencing/run{run}/") for run in runs]

navbar =  html.Div([
                    dbc.NavbarSimple([
                        dbc.NavItem(dbc.NavLink("Global Updates", href="/apps/global_updates")),
                        dbc.NavItem(dbc.NavLink("Kenya Summary",href="/apps/sequencing/statistics/")),
                        
                        dbc.DropdownMenu(dropdown_items, nav=True, in_navbar=True, label="Sequencing Runs")

                        #dbc.NavItem(dbc.NavLink("Submissions", href="/apps/submissions")),
                        #dbc.NavItem(dbc.NavLink("Labs", href = "/apps/labs")),
                    ],              brand_href="/apps/global_updates", 
                                    brand="CIMEA-VEC Reports",
                                    style={"margin-bottom":5},
                                    color="#333972",dark=True,light=True,
                                    fixed ="top",className = "text-light font-weight-bold"
                                    )          
]) 

app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    navbar,
    html.Div(id='page-content', children=[])
])
@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == "/apps/global_updates":
        return global_updates.layout
    elif pathname == "/apps/sequencing/statistics/":
        return sequencing.statistics_layout
    elif pathname == "/apps/sequencing/run231/":
        return run231.updates_layout
    elif pathname == "/apps/sequencing/run232/":
        return run232.updates_layout
    elif pathname == "/apps/sequencing/run233/":
        return run233.updates_layout
    elif pathname == "/apps/sequencing/run234/":
        return run234.updates_layout
    elif pathname == "/apps/sequencing/run235/":
        return run235.updates_layout
    elif pathname == "/apps/sequencing/run236/":
        return dcc.Loading(run236.updates_layout)
    elif pathname == "/apps/sequencing/run237/":
        return dcc.Loading(run237.updates_layout)
    elif pathname == "/apps/sequencing/run238/":
        return dcc.Loading(run238.updates_layout)
    elif pathname == "/apps/sequencing/run239/":
        return dcc.Loading(run239.updates_layout)
    
    elif pathname == "/apps/sequencing/reports/":
        return sequencing.updates_layout
    else:
        return global_updates.layout
        
if __name__ == '__main__':
    app.run_server(debug=True,host="0.0.0.0", port = "8856") #8888