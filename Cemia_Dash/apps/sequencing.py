import dash
import dash_bootstrap_components as dbc
from dash import html, dcc, dash_table
import pandas as pd
import pathlib
import base64
import plotly.express as px
from app import app
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go

#import dataset
PATH = pathlib.Path(__file__).parent
DATA_PATH = PATH.joinpath("../data/").resolve()

markercolor = "#8B0000"
color_patterns = ["#FF5733","#8E44AD","#2236A0","#252525","#1B6311"]
pcolor = "#FFFAFA"
pcolor_white = "white"
gridcolor="lightgray"
pcolor_home = "#E6E6E6"
cardbody_style = {"background-color":pcolor}
margin = dict(l=20, r=25, t=20, b=20)

labs = pd.read_table(DATA_PATH.joinpath("submitting_labs.txt"))
count_ns = pd.read_csv(DATA_PATH.joinpath("count_ns.csv"))
variants_kenya = pd.read_table(DATA_PATH.joinpath("variant_data_kenya.tsv"))
variants_kenya["Month"] = pd.to_datetime(variants_kenya["Month"], format = "%Y-%m-%d")
omicron = pd.read_table(DATA_PATH.joinpath("omicron_data.tsv"))
omicron["Month"] = pd.to_datetime(omicron["Month"], format = "%Y-%m-%d")

#values
variant_values = variants_kenya.groupby('Month')[["Frequency"]].sum().reset_index()
#kwtrp submission
total = len(labs)
kwtrp_no = len(labs[labs["submitting_lab"] == "KEMRI-Wellcome Trust Research Programme,Kilifi"])
percentage_kwtrp = round(kwtrp_no/total*100)
#encoded_image = base64.b64encode(open(image_filename, 'rb').read())
##views Ns

fig_ns = px.scatter(count_ns, x = "date_submitted",y="n_percentage", labels = {"n_percentage":"Proportion of Ns"},range_y = [0,50])
fig_ns.update_layout(height =  400, width=600,plot_bgcolor = pcolor_white)
fig_ns.update_traces(marker_size=2,marker_color = markercolor)
fig_ns.update_xaxes(title= "submission date",linecolor = "black",ticks="outside",nticks=10,tickfont = dict(size=8),title_font = {"size":10})
fig_ns.update_yaxes(linecolor = "black",ticks="outside",tickfont = dict(size=8),title_font = {"size":10},gridcolor = gridcolor)

kwtrp_ns = count_ns.loc[count_ns["submitting_lab"] == "KEMRI-Wellcome Trust Research Programme,Kilifi"]
fig_kwtrp_ns = px.scatter(kwtrp_ns, x = "date_submitted", y = "n_percentage",labels = {"n_percentage":"Proportion of Ns"},range_y = [0,50])
fig_kwtrp_ns.update_layout(height =  400, width=600, plot_bgcolor = pcolor_white)
fig_kwtrp_ns.update_traces(marker_size=2,marker_color = markercolor)
fig_kwtrp_ns.update_xaxes(title= "submission date",linecolor = "black",ticks="outside",tickfont = dict(size=8),title_font = {"size":10}) #dtick="M1",tickformat="%b\n%Y",
fig_kwtrp_ns.update_yaxes(linecolor = "black",ticks="outside",tickfont = dict(size=8),title_font = {"size":10},gridcolor = gridcolor)


#variant plot
fig_var = px.bar(variants_kenya, x = "Month", y = "percentage", color="variant",range_x=["2020-01-01","2023-01-12"],
                 color_discrete_sequence = ["#1b9e77","#d95f02","#7570b3","#e7298a","#8111A5","#e6ab02","#a6761d"])
#fig_var.update_traces(textfont_size=10,textposition="outside", text = variant_values["Frequency"])
fig_var.update_xaxes(nticks = 10,linecolor = "black",ticks="outside",tickfont = dict(size=10),title = None)
fig_var.update_yaxes(linecolor = "black",ticks="outside",tickfont = dict(size=12),title ="Proportion", title_font = {"size":12})
fig_var.update_layout(legend=dict(itemsizing="constant",title = None,orientation = "h", font=dict(size=10)),plot_bgcolor = pcolor_white,
                      margin = margin)


#omicron variant
fig_omicron = px.bar(omicron, x = "Month", y = "percentage", color="pangolin_lineage",
                     color_discrete_sequence = px.colors.qualitative.Prism,range_x=["2021-10-01","2022-12-31"])
fig_omicron.update_yaxes(linecolor = "black",ticks="outside",tickfont = dict(size=10),title ="Proportion", title_font = {"size":12})
fig_omicron.update_xaxes(linecolor = "black",ticks="outside",tickfont = dict(size=10),title =None, title_font = {"size":12})
fig_omicron.update_layout(legend=dict(itemsizing="constant",title = None,orientation = "h", font=dict(size=10)),plot_bgcolor = pcolor_white,
                      margin = margin)

#setting the layout
statistics_layout = html.Div([
    dbc.Spinner([
    dbc.Row([
        html.H5("A summary of sequencing data and process at KWTRP"),
        dbc.Col([
            dbc.CardBody([
                html.H2("12,180"),html.A("Genomes from Kenya")
            ])
        ]),
        dbc.Col([
            dbc.CardBody([
                html.H2(f"{percentage_kwtrp}%"),html.A("submitted by KWTRP")
            ])
        ]),
    ],justify="evenly",className = "mt-5 pt-5 ms-5 ps-4"),
    
    html.Hr(className = "ms-4 me-4"),
    
    dbc.Row([
        
        dbc.Col([
            html.Label("Temporal distribution of SARS-CoV-2 variants in Kenya",style = {"text-align":"center","font-size":14}, className = "text-center fw-bold text-dark"),
            dbc.Spinner([dcc.Graph(figure = fig_var, responsive = True,style = {"height":"300px", "width":"800px"})]),
        ],width = {"size":8,"offset":1},lg=9) 
           
    ],justify="center"),
    
    html.Hr(className = "ms-4 me-4"),
    dbc.Row([
        
        dbc.Col([
            html.Label("SARS-CoV-2 Omicron lineages in Kenya (November 2021 - December 2022)",style = {"text-align":"center","font-size":14}, className = "fs-6 text-center text-dark"),
            dcc.Graph(figure = fig_omicron, responsive = True,style = {"width":"700px","height":"300px"})
        ],width = {"size":8,"offset":1},lg=9)    
    ],justify="center"),
    
    
    
    
    html.Hr(className = "ms-4 me-4"),
    
    dbc.Row([
        dbc.Col([
            html.P("Proportion of Ns in sequences across Kenya",className = "text-center fw-bold fs-6 text-dark"), #fs-6 fst-italic
            dcc.Graph(figure = fig_ns),
        ],width=5),
        dbc.Col([
            html.P("Proportion of Ns in sequences from KWTRP",className = "text-center fw-bold fs-6 text-dark"), 
            dcc.Graph(figure=fig_kwtrp_ns)
        ],width=5)
    
    ]),
    
    html.Hr(className = "ms-4 me-2"),
    ])
    
    # dbc.Row([
    #     dbc.Col([
    #         html.H6("Average coverage (breadth) of sequences submitted to GISAID database."),
    #         html.Img(src = app.get_asset_url("coverage_plot.png"))
    #     ],width = 5),
    #     dbc.Col([
            
    #     ],width=5)
    # ],className = "ps-4 pe-4")
])

#