import dash
import dash_bootstrap_components as dbc
from dash import html, dcc, dash_table
import pandas as pd
import pathlib
from dash_bootstrap_templates import load_figure_template
import base64
import plotly.express as px
from app import app
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from utils import *

#load_figure_template("pulse")
#import dataset
PATH = pathlib.Path(__file__).parent
DATA_PATH = PATH.joinpath("../data").resolve()


col_title = "text-center text-black fw-normal"
markercolor = "#8B0000"
margin = dict(l=10, r=10, t=5, b=5)
color_patterns = ["#FF5733","#8E44AD","#2236A0","#252525","#1B6311"]
pcolor = "rgba(0,0,0,0)"
pcolor_white = "white"
gridcolor="lightgray"
pcolor_home = "#E6E6E6"
cardbody_style = {"background-color":pcolor}
style_text ={"font-size":14,"text-align":"center"}

complete_data = pd.read_table(DATA_PATH.joinpath("run234/complete_data.txt"),index_col = "sample_id")
sample_source = pd.read_table(DATA_PATH.joinpath("run234/sample_source.txt"))

tickfont = dict(size=8)
titlefont = {"size":10}

total_samples = len(complete_data)
samples_sequenced = len(complete_data[complete_data["sequence"].notna()])
met_threshold =  len(complete_data[(complete_data["sequence"].notna()) & (complete_data["coverage"] >= 0.7)])
min_ns = int(round(complete_data["n_proportion"].min() * 100,0))
max_ns = int(round(complete_data["n_proportion"].max() * 100,0))

###Box plot for cycle threshold distribution view. 
# fig_box = px.box(pre_data, y = ["failed_pcr","not_sequenced","failed_sequencing"],points="all",
#                  color_discrete_sequence = ["#8B0000"],range_y = [0,40])
# fig_box.update_traces(width=0.1, marker_size=4)
# fig_box.update_layout(uniformtext_minsize=10, margin = margin,paper_bgcolor = pcolor, plot_bgcolor = pcolor)
# fig_box.update_yaxes(title =  "Cycle Threshold (ct)",linecolor = "black",ticks="outside",tickfont = tickfont,
#                      title_font = titlefont,
#                      gridcolor="lightgray",linewidth=1,nticks=10)
# fig_box.update_xaxes(linecolor = "black",ticks="outside", tickfont = tickfont,title=None)


#sequence coverage by sequence length
cov_data = complete_data[complete_data["coverage"].notna()].sort_values("coverage",ascending=True)
fig_cov = px.bar(cov_data, x = cov_data.index, y = "coverage")
fig_cov.update_layout(uniformtext_minsize=8,plot_bgcolor = pcolor,paper_bgcolor = pcolor,margin=margin)
fig_cov.update_xaxes(title =None,linecolor = "black", tickfont = tickfont,showticklabels = False,title_font = titlefont)
fig_cov.update_yaxes(title = "Coverage",linecolor = "black", tickfont = tickfont,title_font = titlefont)
fig_cov.update_traces(marker_color = "#51A2C4")
fig_cov.add_hrect(y0=0.70, y1=0.70, line_width=1, fillcolor="red", opacity=0.50)

##Proportion of Ns
#prop_recovered = 100 - round(max(lineage["n_proportion"])*100,0)

fig_n_prop = px.scatter(complete_data.sort_values("n_proportion",ascending=False), range_y = [0,1],x = complete_data.index, 
                        y = "n_proportion",labels = {"n_proportion":"Proportion of Ns (%)"})
fig_n_prop.update_yaxes(title_font = titlefont,linecolor = "black",tickfont = tickfont)
fig_n_prop.update_xaxes(title = None,title_font = titlefont,linecolor = "black",showticklabels = False)
fig_n_prop.update_traces(marker_color = "#51A2C4",marker_size=4)

# fig_n_prop.update_yaxes(linecolor = "black", tickfont = tickfont,title_font = titlefont,gridcolor = gridcolor)
# fig_n_prop.update_xaxes(title = "sample",linecolor = "black", tickfont = tickfont,showticklabels = False,
#                         title_font = titlefont)
fig_n_prop.update_layout(margin = margin,plot_bgcolor = pcolor,paper_bgcolor = pcolor)


#observed variants
unassigned = len(complete_data[complete_data["lineage"] == "Unassigned"])
variants_data = (complete_data[(complete_data["sequence"].notna()) & (complete_data["lineage"] != "Unassigned")])
variants_data = variants_data.groupby("lineage")[["lineage"]].count().rename(columns = {"lineage":"Freq"}).sort_values("Freq",ascending=False)
variants_data.sort_values("Freq", ascending=True,inplace=True)
fig_var = px.bar(variants_data, x = "Freq",y = variants_data.index,text_auto=True)
fig_var.update_traces(width = 0.6,textfont_size=8,textposition="outside",marker_color = "#51A2C4",cliponaxis=True)
fig_var.update_layout(uniformtext_minsize=8,margin=margin,bargap=0.08,plot_bgcolor = pcolor,paper_bgcolor = pcolor)
fig_var.update_yaxes(linecolor = "black",ticks="outside", tickfont = tickfont,title=None)
fig_var.update_xaxes(title = None,linecolor = "black",ticks="outside", tickfont = {"size":8},title_font = titlefont)



#coverage vs ct-value
ct_cov = complete_data.loc[:,["coverage", "ctvalues"]]
ct_cov.sort_values("ctvalues", ascending=True, inplace = True)
fig_ct_cov = px.scatter(ct_cov, y = "coverage",x = "ctvalues", trendline = "ols", range_y = [0,1],range_x = [0,35])#,trendline_options=dict(log_x=True))
fig_ct_cov.update_yaxes(title = "Coverage (breadth)",title_font = titlefont, tickfont = tickfont, linecolor="black")
fig_ct_cov.update_xaxes(gridcolor = gridcolor,nticks = 10,title = "CT-Value (E gene)",title_font = titlefont, 
                        tickfont = tickfont, linecolor="black")
fig_ct_cov.update_traces(marker_color = "#51A2C4",marker_size=4)
fig_ct_cov.update_layout(plot_bgcolor = pcolor,paper_bgcolor = pcolor,margin=margin)



#bubble plot
bubble_data =  complete_data.loc[:,["date_collected","lineage","Facility"]]
bubble_data= pd.DataFrame(bubble_data.groupby(["date_collected","Facility","lineage"])["lineage"].count())\
                .rename(columns = {"lineage":"Freq"}).reset_index()
bubble_data =  bubble_data[bubble_data["lineage"] != "Unassigned"]
bubble_data["date_collected"] = pd.to_datetime(bubble_data["date_collected"]) #, format="%d/%m/%Y"
bubble_data.sort_values("date_collected", inplace = True)

fig_bubble = px.scatter(bubble_data, y = "Facility", x = "date_collected", size_max=15,
                                size = "Freq", color="lineage")
fig_bubble.update_layout(plot_bgcolor = pcolor,margin=margin,paper_bgcolor = pcolor,
                         legend=dict(itemsizing="constant",title = None,borderwidth=0,orientation = "h", font={"size":12},
                                     itemwidth=30, yanchor = "top",y = 1.2,xanchor = "left"))
fig_bubble.update_xaxes(gridcolor = gridcolor,title = None, title_font = dict(size=12),
                        tickfont = tickfont,linecolor = "black",ticks="outside")
fig_bubble.update_yaxes(gridcolor = gridcolor,title = None, title_font = dict(size=12),
                        tickfont = tickfont,linecolor = "black",ticks="outside")


card_style = "bg-light border shadow"# shadow"
classname_col = "bg-light bg-opacity-20 g-1 justify-content-center p-2 m-2" 
style = {"height":"300px", "width":"400px"}


updates_layout = html.Div([
            dbc.Row([
                dbc.Row([
                    dbc.Col([
                        html.Br(),
                        html.H4("Summary of Run234", className='text-center'),
                        html.P([f"In Run 234, a total of ", html.B(total_samples), " samples from ",html.B("Kilifi"),\
                            " were processed.",html.B(samples_sequenced), html.B(f" ({round(samples_sequenced/total_samples*100)}%)") ,\
                                " were successfully sequenced. Of the sequenced, ", html.B(met_threshold),html.B(f" ({round(met_threshold/samples_sequenced*100)}%)" ) , \
                                " met the coverage threshold of >70%", ". \
                                The proportion of Ns ranged between ", html.B(f"{min_ns} - {max_ns}%")],
                        className = "text-start",style={"font-size":15,}),

                    ],width=8,className="bg-light"),
                ], justify = "center",className = "mt-5 pt-5"),

                dbc.Row([

                    dbc.Col([
                        html.Br(),
                        html.P("Sample Summary",className = col_title),
                        dash_table.DataTable(sample_source.to_dict("records"), [{"name":i, "id":i} for i in sample_source.columns],
                                     page_size = 4,
                                     #style_header = {"backgroundColor":"gray",'color': 'white','fontWeight': 'bold',"font-size":12, "text-align":"center"},
                                     #style_data = {"backgroundColor":"white","color":"black","font-size":12},
                                     #style_cell = {"textAlign":"left","minWidth":"5px","maxWidth":"180px"},
                                     style_cell = {"whiteSpace":"normal","height":"auto"},
                                     #style_table = {"height":"900px","overflowX":"auto"},
                                     #editable = True
                                     style_data = {"font-size":10,"text-align":"center"},#"width":"20px", "height":"10px"},
                                     style_header = {"font-size":14,'fontWeight': 'bold',"text-align":"center",
                                                     "backgroundColor":"gray",'color': 'white'},
                                     style_table={"height":"150px","overflowX":"auto","backgroundColor":pcolor}
                                    ),
                    ],width =5,className = "h-100 bg-light  border-bottom", style = {"margin-right":"10px"}),

                    dbc.Col([
                        html.Br(),
                        html.P("Correlation for CT vs Sequence Coverage",className = col_title),
                        dcc.Graph(figure = fig_ct_cov, responsive = True, style = {"height":"30vh","width":"20hw"})
                    ],width=3,className = "bg-light border-bottom"),

                ], justify = "center",align="center", className = "mt-2"),

                dbc.Row([
                    dbc.Col([
                        html.Br(),
                            html.P("Sequence Coverage",className = col_title),
                            dcc.Graph(figure = fig_cov, responsive = True, style = {"height":"30vh", "width":"20hw"}), #"width":"20hw",
                            html.Br(),
                            html.P([f"Sequences >70% coverage get submitted to GISAID. For this run, only ",html.B(met_threshold), " sequences passed the QC \
                                test."],style = style_text),
                            
                    ], width=4,className = "bg-light",style = {"margin-right":"10px"} ), #

                    dbc.Col([ 
                             html.Br(),
                            html.P("Proportion of Ns in Sequences",className = col_title),
                               dcc.Graph(figure = fig_n_prop, responsive = True, style = {"height":"30vh","width":"20hw"}),
                            html.Br(),
                            html.P([f"Genome recovery ranged between ", html.B(f"{min_ns} - {max_ns}%."),  "of the sequence."] ,style = style_text),
                    ],width=4,className = "h-100 bg-light"),        

                ],justify = "center",align="center", className = "mt-2"),

                dbc.Row([
                    dbc.Col([
                        html.Br(),
                        html.P("Variant Frequency",className = col_title),
                        dcc.Graph(figure = fig_var, responsive = True, style = {"height":"25vh"}), #"width":"32hw",
                        html.Br(),
                        html.P([f"Unassigned sequences - ", html.B(unassigned) ], style = style_text),
                    ],width=3,className = "h-100 bg-light",style = {"margin-right":"10px"}),

                    dbc.Col([
                        html.Br(),
                        html.P("Temporal Variant Distribution",className = col_title),
                        dcc.Graph(figure = fig_bubble, responsive = True,style = {"height":"35vh"}), #"width":"32hw"
                        html.Br(),
                        html.P([f"There is countinued detection of ", html.B("XBB.1.22.1"), " assigned as ", html.B("FY.4. "),
                                
                        ],style = style_text),
                    ],width=5,className = "h-100 bg-light"),

                ],justify="center",align="center",className = "mt-2"),
                
                dbc.Row([
                    html.P([f"Lineage classification was performed on ",html.B("May 5, 2023 "),"using", \
                        html.B(" pangolin v.4.2") ], className = "text-dark",style = style_text)
                ],justify="center",align="center",className = "mt-2 mb-4")
                
            ], className = "bg-secondary bg-opacity-25 g-1 justify-content-center ps-4 pe-4 m-2")
             
])
                   