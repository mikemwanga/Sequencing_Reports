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

load_figure_template("flatly")
#import dataset
PATH = pathlib.Path(__file__).parent
DATA_PATH = PATH.joinpath("../data").resolve()

labs = pd.read_table(DATA_PATH.joinpath("submitting_labs.txt"))
markercolor = "#8B0000"
margin = dict(l=20, r=25, t=20, b=20)
color_patterns = ["#FF5733","#8E44AD","#2236A0","#252525","#1B6311"]
pcolor = "rgba(0,0,0,0)"
pcolor_white = "white"
gridcolor="lightgray"
pcolor_home = "#E6E6E6"
cardbody_style = {"background-color":pcolor}

count_ns = pd.read_csv(DATA_PATH.joinpath("count_ns.csv"))
sample_source = pd.read_table(DATA_PATH.joinpath("run231/sample_source.tsv"))
pre_data = pd.read_table(DATA_PATH.joinpath("run231/lab_sample_summary.tsv"),index_col = "sample_id")
lineage = pd.read_table(DATA_PATH.joinpath("run231/lineage_coverage.tsv"),index_col="sample_id")
#lineage.drop("S11649",inplace=True)
lineage.sort_values("coverage", ascending = False,inplace=True)


tickfont = dict(size=12)
titlefont = {"size":14}

failed_submission = len(lineage[lineage["coverage"] < 0.7])
total_samples = sample_source["#samples"].sum()
site = sample_source["Facility"]
samples_sequenced = sample_source["sequenced"]
failed_prop = round(failed_submission/len(lineage)*100,1)

###Box plot for cycle threshold distribution view. 

fig_box = px.box(pre_data, y = ["failed_pcr","not_sequenced","failed_sequencing"],points="all",
                 color_discrete_sequence = ["#8B0000"],range_y = [0,40])
fig_box.update_traces(width=0.1, marker_size=4)
fig_box.update_layout(uniformtext_minsize=10, margin = margin,paper_bgcolor = pcolor, plot_bgcolor = pcolor)
fig_box.update_yaxes(title =  "Cycle Threshold (ct)",linecolor = "black",ticks="outside",tickfont = tickfont,
                     title_font = titlefont,
                     gridcolor="lightgray",linewidth=1,nticks=10)
fig_box.update_xaxes(linecolor = "black",ticks="outside", tickfont = tickfont,title=None)


#sequence coverage by sequence length
fig_cov = px.bar(lineage, x = lineage.index, y = "coverage")
fig_cov.update_layout(uniformtext_minsize=8,plot_bgcolor = pcolor,paper_bgcolor = pcolor,margin=margin)
fig_cov.update_xaxes(title = "sample",linecolor = "black", tickfont = tickfont,showticklabels = False,title_font = titlefont)
fig_cov.update_yaxes(title = "Coverage",linecolor = "black", tickfont = tickfont,title_font = titlefont)
fig_cov.update_traces(marker_color = "#51A2C4")
fig_cov.add_hrect(y0=0.70, y1=0.70, line_width=1, fillcolor="red", opacity=0.50)

##Proportion of Ns
prop_recovered = 100 - round(max(lineage["n_proportion"])*100,0)

fig_n_prop = px.scatter(lineage.sort_values("n_proportion",ascending=False), range_y = [0,1],x = lineage.index, 
                        y = "n_proportion",labels = {"n_proportion":"Proportion of Ns (%)"})
fig_n_prop.update_traces(marker_color = "#51A2C4",marker_size=4)
fig_n_prop.update_yaxes(linecolor = "black", tickfont = tickfont,title_font = titlefont,gridcolor = gridcolor)
fig_n_prop.update_xaxes(title = "sample",linecolor = "black", tickfont = tickfont,showticklabels = False,
                        title_font = titlefont)
fig_n_prop.update_layout(plot_bgcolor = pcolor,paper_bgcolor = pcolor,margin = margin)

#observed variants
variants = lineage.groupby("Nextclade_pango")[["Nextclade_pango"]].count().rename(columns = {"Nextclade_pango":"Freq"}).sort_values("Freq",ascending=False)
variants.sort_values("Freq", ascending=True,inplace=True)
fig_var = px.bar(variants, x = "Freq",y = variants.index,text_auto=True)
fig_var.update_traces(width = 0.6,textfont_size=8,textposition="outside",marker_color = "#51A2C4",cliponaxis=True)
fig_var.update_layout(uniformtext_minsize=8,margin=margin,bargap=0.08,plot_bgcolor = pcolor,paper_bgcolor = pcolor)
fig_var.update_yaxes(linecolor = "black",ticks="outside", tickfont = tickfont,title=None)
fig_var.update_xaxes(title = "Frequency",linecolor = "black",ticks="outside", tickfont = {"size":8},title_font = titlefont)

#coverage vs ct-value
ct_cov = lineage.loc[:,["coverage", "E MetaBion"]]
ct_cov.sort_values("E MetaBion", ascending=True, inplace = True)
fig_ct_cov = px.scatter(ct_cov, y = "coverage",x = "E MetaBion", trendline = "ols", range_y = [0,1],range_x = [0,35])#,trendline_options=dict(log_x=True))
fig_ct_cov.update_yaxes(title = "Coverage (breadth)",title_font = titlefont, tickfont = tickfont, linecolor="black")
fig_ct_cov.update_xaxes(gridcolor = gridcolor,nticks = 10,title = "CT-Value (E gene)",title_font = titlefont, 
                        tickfont = tickfont, linecolor="black")
fig_ct_cov.update_traces(marker_color = "#51A2C4",marker_size=4)
fig_ct_cov.update_layout(plot_bgcolor = pcolor,paper_bgcolor = pcolor,margin=margin)

#site based comparisons
facility_lineage = lineage.loc[:,["Nextclade_pango","facility"]].groupby(["Nextclade_pango","facility"])[["Nextclade_pango"]].count().rename(columns = {"Nextclade_pango":"Freq"}).reset_index()
unassigned = facility_lineage[facility_lineage["Nextclade_pango"] != "Unassigned"]
facility_lineage =  facility_lineage[facility_lineage["Nextclade_pango"] != "Unassigned"]
no_facilities = len(facility_lineage["facility"].unique())
plot_title  = facility_lineage["facility"].unique()
fig_site = make_subplots(rows = no_facilities, cols=1,shared_xaxes=True, horizontal_spacing = 0.8, subplot_titles=plot_title)
row = 1
for f in facility_lineage["facility"].unique():
    data = facility_lineage[facility_lineage["facility"] == f]
    fig_site.add_trace(go.Bar(x=data["Nextclade_pango"], y = data["Freq"], name = f,showlegend=False),
                row=row, col=1)
    row += 1
fig_site.update_yaxes(range = [0,30],tickfont = tickfont)
fig_site.update_xaxes(tickfont = dict(size=8), tickangle = -60)
fig_site.update_traces(textfont_size=8,width = 0.4)
fig_site.update_layout(font=dict(size=10,color="black"),margin=margin)#,title_font=tickfont)#,title = dict(font=dict(size=8)))
for i in fig_site["layout"]["annotations"]:
    i["font"] = tickfont

#bubble plot
bubble_data = lineage.loc[:,["date_collected","Nextclade_pango","facility"]]
bubble_data= pd.DataFrame(bubble_data.groupby(["date_collected","facility","Nextclade_pango"])["Nextclade_pango"].count()).rename(columns = {"Nextclade_pango":"Freq"}).reset_index()
bubble_data =  bubble_data[bubble_data["Nextclade_pango"] != "Unassigned"]
bubble_data["date_collected"] = pd.to_datetime(bubble_data["date_collected"], format="%d/%m/%Y")
bubble_data.sort_values("date_collected", inplace = True)
fig_bubble = px.scatter(bubble_data, y = "facility", x = "date_collected", size_max=30,size = "Freq", color="Nextclade_pango")
fig_bubble.update_layout(plot_bgcolor = pcolor,margin=margin,paper_bgcolor = pcolor,
                         legend=dict(itemsizing="constant",title = None,borderwidth=0,orientation = "h", font={"size":12},
                                     itemwidth=30, yanchor = "top",y = 1.2,xanchor = "right",x = 1))
fig_bubble.update_xaxes(gridcolor = gridcolor,title = None, title_font = dict(size=12),
                        tickfont = tickfont,linecolor = "black",ticks="outside")
fig_bubble.update_yaxes(gridcolor = gridcolor,title = None, title_font = dict(size=12),
                        tickfont = tickfont,linecolor = "black",ticks="outside")

card_style = "bg-light border shadow"# shadow"
classname_col = "bg-light bg-opacity-20 g-1 justify-content-center p-2 m-2" 
style = {"height":"300px", "width":"400px"}

updates_layout = html.Div([
    dbc.Row([
        html.H4("A summary of processed lab samples for Run 231",className = "text-center fs-5 text-dark"),
        dbc.Row([
        dbc.Col([
        html.P(["In Run 231, a total of ", html.B(total_samples), " samples from ", html.B(site), " were processed and ",html.B(samples_sequenced), \
            " were successfully sequenced and processed. Only ", html.B("2"), " \
            met the coverage threshold of ", html.B(">70%"),". The proportion of Ns ranged \
            between ", html.B("0-28%"), ". The following lineages were identified: ", html.B("XBB.1.5(n=2).")],
            className = "text-start",
            style={"font-size":15,}),
        ], width = 9,lg=8),
        ],justify = "center"),
        
        dbc.Col([
            html.P("Sources of samples processed",style = {"text-align":"center","font-size":titlefont},
                                className = "text-center fw-bold text-dark"),
            html.Br(),html.Br(),
            dash_table.DataTable(sample_source.to_dict("records"), [{"name":i, "id":i} for i in sample_source.columns],
                                 page_size = 4,
                                 #style_header = {"backgroundColor":"gray",'color': 'white','fontWeight': 'bold',"font-size":12, "text-align":"center"},
                                 #style_data = {"backgroundColor":"white","color":"black","font-size":12},
                                 #style_cell = {"textAlign":"left","minWidth":"5px","maxWidth":"180px"},
                                 style_cell = {"whiteSpace":"normal","height":"auto"},
                                 #style_table = {"height":"900px","overflowX":"auto"},
                                 #editable = True
                                 style_data = {"font-size":13,"text-align":"center"},#"width":"20px", "height":"10px"},
                                 style_header = {"font-size":14,'fontWeight': 'bold',"text-align":"center",
                                                 "backgroundColor":"gray",'color': 'white'},
                                 style_table={"height":"150px","overflowX":"auto","backgroundColor":pcolor}
                                ),

            html.Br(),html.Br(),html.Br(),html.Br(), html.Hr(),
            
            #another graph
            html.P("Coverage for sequenced samples",style = {"text-align":"center","font-size":titlefont},
                       className = "text-center fw-bold text-dark"),
            dcc.Graph(figure = fig_cov, responsive = True, style = {"height":"300px", "width":"300px"},),
            
            html.Label("Sequences >70% coverage get submitted to GISAID database. For this run, all sequenceses passed the QC test.",
                       style = {"text-align":"start","font-size":14}),
            html.Label(f"{failed_prop}%({failed_submission}) failed submission QC",style = {"text-align":"center","font-size":14}),
            
            html.Hr(),
            html.P("Correlation of CT values to sequence coverage (breadth)",
                    style = {"text-align":"center","font-size":titlefont},className = "fw-bold text-dark"),
             dcc.Graph(figure = fig_ct_cov, responsive = True, style = {"height":"300px", "width":"350px"})
            
            
        ],width=5,xl=4,className = card_style, align = "center",style = {"height":"150vh"}),
        
        dbc.Col([
            html.P("Summary of samples that failed PCR and pre-sequencing process",className = "text-center fw-bold text-dark",
                       style = {"text-align":"center","font-size":titlefont}),
            
            dcc.Graph(figure = fig_box, responsive = True, style = {"height":"280px", "width":"350px"},className = "ms-3"),
            html.Hr(),
            
            #graph2
            html.P("Proportion of Ns in sequences",
                    style = {"text-align":"center","font-size":titlefont},className = "fw-bold text-dark"),
            dcc.Graph(figure = fig_n_prop, responsive = True, style = {"height":"300px", "width":"350px"},className = "ms-3"),

            html.Label(f"Genome recovery was over {prop_recovered}%",
                       style = {"text-align":"center","font-size":14}),
            html.Br(),html.Hr(),
            
            html.P("SARS-CoV-2 variants from the current run",
                    style = {"text-align":"center","font-size":titlefont},className = "fw-bold text-dark"),
            dcc.Graph(figure = fig_var, responsive = True, style = {"height":"300px", "width":"400px"}),
            html.Label("BQ.1.1 is the most abundant lineage at 73%",
                       style = {"text-align":"center","font-size":14}),
        ],width=5,xl=4,className = card_style, style = {"height":"150vh", "margin-left":"5px"})
        
        
    ],justify="center",className = "mt-5 pt-5 bg-secondary bg-opacity-10",align = "center"),
    
    dbc.Row([
        html.Br(),html.Br(),
        dbc.Col([
            html.Br(),
            html.P("Temporal visualization of isolated SARS-CoV-2 variants",
                    style = {"font-size":titlefont},className = "text-center fw-bold text-dark"),
            dcc.Graph(figure = fig_bubble, responsive = True,style = {"height":"350px", "width":"800px"})
        ],width=9,lg=8, className = card_style),
        
        
    ],justify="center",className = "bg-secondary bg-opacity-10 mt-3 mb-3"),
    
])