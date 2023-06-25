
from utils import *

#load_figure_template("pulse")
#import dataset
PATH = pathlib.Path(__file__).parent
DATA_PATH = PATH.joinpath("../data").resolve()


runid="run236"

col_title = "text-center text-black fw-normal"
markercolor = "#8B0000"
margin = dict(l=10, r=10, t=5, b=5) #'#543005'
color_patterns = ['#1b9e77','#d95f02','#7570b3','magenta','#3B0B0B','#40004b','#543005','#053061','#006837','#6a3d9a',"#67001f"]
pcolor = "rgba(0,0,0,0)"
pcolor_white = "white"
gridcolor="lightgray"
pcolor_home = "#E6E6E6"
cardbody_style = {"background-color":pcolor}
style_text ={"font-size":14,"text-align":"center"}

complete_data = pd.read_table(DATA_PATH.joinpath(runid+"/complete_data.txt"),index_col = "sample_id")
sample_source = pd.read_table(DATA_PATH.joinpath(runid+"/sample_source.txt"))

variants_data_2 = pd.read_table(DATA_PATH.joinpath("kenya_lineages.txt"))
variant_map = pd.read_table(DATA_PATH.joinpath("map_lineages.tsv"))

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
unassigned = len(complete_data[complete_data["Nextclade_pango"] == "Unassigned"])
variants_data = (complete_data[(complete_data["sequence"].notna()) & (complete_data["Nextclade_pango"] != "Unassigned") & (complete_data["coverage"] >=0.7)])
variants_data_group = variants_data.groupby("Nextclade_pango")[["Nextclade_pango"]].count().rename(columns = {"Nextclade_pango":"Freq"}).sort_values("Freq",ascending=False)
variants_data_group.sort_values("Freq", ascending=True,inplace=True)
fig_var = px.bar(variants_data_group, x = "Freq",y = variants_data_group.index,text_auto=True)
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
bubble_data =  variants_data.loc[:,["date_collected","Nextclade_pango","Facility"]]
bubble_data= pd.DataFrame(bubble_data.groupby(["date_collected","Facility","Nextclade_pango"])["Nextclade_pango"].count())\
                .rename(columns = {"Nextclade_pango":"Freq"}).reset_index()
bubble_data =  bubble_data[bubble_data["Nextclade_pango"] != "Unassigned"]
bubble_data["date_collected"] = pd.to_datetime(bubble_data["date_collected"],format="%d/%m/%Y") #, format="%d/%m/%Y"
bubble_data.sort_values("date_collected", inplace = True)

#fig_bubble = px.bar(bubble_data, x = "Facility", y = "Freq"
#                                , color="Nextclade_pango",color_discrete_sequence=color_patterns)
fig_bubble = px.histogram(bubble_data, x = "Facility", y = "Freq"
                                , color="Nextclade_pango",barmode="group",color_discrete_sequence=color_patterns)
#fig_bubble = px.scatter(bubble_data, y = "Facility", x = "date_collected", size_max=15,
#                                size = "Freq", color="Nextclade_pango",color_discrete_sequence=color_patterns)
fig_bubble.update_layout(plot_bgcolor = pcolor,margin=margin,paper_bgcolor = pcolor,
                         legend=dict(itemsizing="constant",title = None,borderwidth=0,orientation = "h", font={"size":10},
                                     itemwidth=30, yanchor = "bottom",y = -0.5,xanchor = "left"))
fig_bubble.update_xaxes(gridcolor = None,title = None, title_font = dict(size=14),
                        tickfont = dict(size=12),linecolor = "black",ticks="outside")
fig_bubble.update_yaxes(gridcolor = gridcolor,title = None, title_font = dict(size=12),
                        tickfont = tickfont,linecolor = "black",ticks="outside")


lineages = variants_data_2[["date","pangolin_lineage"]].set_index("date")#.head()
lineages.index = pd.to_datetime(lineages.index, format='%d/%m/%Y')
recent_lineages = lineages[lineages.index >= "2023-01-01"].reset_index()
#lineage_map = dict(zip(variant_map.lineage, variant_map.lineage_group))
#recent_lineages["lineage_group"] = recent_lineages["pangolin_lineage"].map(lineage_map)
recent_lineages = recent_lineages[recent_lineages["pangolin_lineage"] != "Unassigned"] #drop columns with unassigned annotations
recent_lineages = recent_lineages.groupby(["date","pangolin_lineage"])[["pangolin_lineage"]].count().\
    rename(columns = {"pangolin_lineage":"Frequency"}).reset_index()

#do the plot
sars_lineages = px.scatter(recent_lineages, x = "date",y = "pangolin_lineage",size= "Frequency",color = "pangolin_lineage",
                    color_discrete_sequence= ['#67001f'],#color_patterns, #
                    range_x=["2023-01-01","2023-06-01"])
sars_lineages.update_layout(plot_bgcolor = pcolor,margin=margin,paper_bgcolor = pcolor,showlegend = False)
sars_lineages.update_xaxes(title = None,linecolor = "black",tickfont = dict(size=10), nticks=6)
sars_lineages.update_yaxes(title = None, linecolor = "black",tickfont = dict(size=10),gridcolor = gridcolor)


card_style = "bg-light border shadow"# shadow"
classname_col = "bg-light bg-opacity-20 g-1 justify-content-center p-2 m-2" 
style = {"height":"300px", "width":"400px"}


updates_layout = html.Div([
            dbc.Row([
                dbc.Row([
                    dbc.Col([
                        html.Br(),
                        html.H4(f"Summary of {runid}", className='text-center'),
                        html.P([f"In {runid}, a total of ", html.B(total_samples), " samples from ",html.B("Kilifi, Lancet and Kwale "),\
                            " were processed.",html.B(f" {samples_sequenced}"), html.B(f" ({round(samples_sequenced/total_samples*100)}%)") ,\
                                " were successfully sequenced. Of the sequenced, ", html.B(met_threshold),html.B(f" ({round(met_threshold/samples_sequenced*100)}%)" ) , \
                                " met the coverage threshold of >70%", ". \
                                The proportion of Ns ranged between ", html.B(f"{min_ns} - {max_ns}%")],
                        className = "text-start",style={"font-size":15,}),

                    ],width=10,className="bg-light"),
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
                                     style_data = {"font-size":12,"text-align":"center"},#"width":"20px", "height":"10px"},
                                     style_header = {"font-size":14,'fontWeight': 'bold',"text-align":"center",
                                                     "backgroundColor":"gray",'color': 'white'},
                                     style_table={"height":"150px","overflowX":"auto","backgroundColor":pcolor}
                                    ),
                    ],width =5,className = "h-100 bg-light  border-bottom", style = {"margin-right":"10px"}),

                    dbc.Col([
                        html.Br(),
                        html.P("Correlation for CT vs Sequence Coverage",className = col_title),
                        dcc.Graph(figure = fig_ct_cov, responsive = True, style = {"height":"30vh","width":"20hw"})
                    ],width=5,className = "bg-light border-bottom"),

                ], justify = "center",align="center", className = "mt-2"),

                dbc.Row([
                    dbc.Col([
                        html.Br(),
                            html.P("Sequence Coverage",className = col_title),
                            dcc.Graph(figure = fig_cov, responsive = True, style = {"height":"30vh", "width":"20hw"}), #"width":"20hw",
                            html.Br(),
                            html.P([f"Sequences >70% coverage get submitted to GISAID. For this run, only ",html.B(met_threshold), " sequences passed the QC \
                                test."],style = style_text),
                            
                    ], width=5,className = "h-100 bg-light",style = {"margin-right":"10px"} ), #

                    dbc.Col([ #fig_n_prop
                             html.Br(),
                            html.P("Variant Frequency",className = col_title),
                               dcc.Graph(figure = fig_var, responsive = True, style = {"height":"30vh","width":"20hw"}),
                            html.Br(),
                            html.P([f"Unassigned sequences - ", html.B(unassigned) ], style = style_text),
                            html.P([f"There is countinued dominance of ", html.B("FY.4. "), "First detection of ",\
                                html.B("XBB.1.22.2, XBB.1.9.2, XBB.1.5.7, and XBB.1.16"),                             
                            ],style = style_text),
                            #html.P([f"Genome recovery ranged between ", html.B(f"{min_ns} - {max_ns}%."),  "of the sequence."] ,style = style_text),

                            #html.P([f"Genome recovery ranged between ", html.B(f"{min_ns} - {max_ns}%."),  "of the sequence."] ,style = style_text),
                    ],width=5,className = "h-100 bg-light"),        

                ],justify = "center",align="center", className = "mt-2"),

                dbc.Row([
                    # dbc.Col([
                    #     html.Br(),
                    #     html.P("Variant Frequency",className = col_title),
                    #     dcc.Graph(figure = fig_var, responsive = True, style = {"height":"25vh"}), #"width":"32hw",
                    #     html.Br(),
                    #     html.P([f"Unassigned sequences - ", html.B(unassigned) ], style = style_text),
                    #     html.P([f"There is countinued detection of ", html.B("FY.4. "), "First detection of ",html.B("XBB.1.22.2,XBB.1.9.2,XBB.1.5.7,XBB.1.16, \
                    #         ,BA.2.23 and EG.1. "), "EG.1 was seen first in USA\
                    #         in January 2023."                             
                    #         ],style = style_text),
                    # ],width=3,className = "h-100 bg-light",style = {"margin-right":"10px"}),
                    
                    dbc.Col([
                        html.Br(),
                        html.P("Variants by Site",className = col_title),
                        dcc.Graph(figure = fig_bubble, responsive = True,style = {"height":"35vh"}), #"width":"32hw"
                        html.Br(),
                        #html.P([f"There is countinued detection of ", html.B("XBB.1.22.1"), " assigned as ", html.B("FY.4. ")],style = style_text),
                    ],width=3,className = "h-100 bg-light"),

                    dbc.Col([
                        html.Br(),
                        html.P("Variant Distribution in Kenya: Jan - May 2023",className = col_title),
                        dcc.Graph(figure = sars_lineages, responsive = True,style = {"height":"35vh"}), #"width":"32hw"
                        html.Br(),
                        
                    ],width=7,className = "h-100 bg-light",style = {"margin-left":"10px"}),

                ],justify="center",align="center",className = "mt-2"),
                
                # dbc.Row([
                #    dbc.Col([
                #        html.P("SARS-COV-2 lineages between January to April 2023",className = col_title),
                #        dcc.Graph(figure = sars_lineages,responsive = True,style = {"width":"40hw","height":"30vh"}),
                #    ],width=8,className = "h-100 bg-light") 
                # ],justify="center",align="center",className = "mt-2"),
                
                dbc.Row([
                    html.P([f"Lineage classification was performed on ",html.B("May 18, 2023 "),"using", \
                        html.B(" pangolin v.4.3") ], className = "text-dark",style = style_text)
                ],justify="center",align="center",className = "mt-2 mb-4")
                
            ], className = "bg-secondary bg-opacity-10 g-1 justify-content-center ps-4 pe-4 m-2")
             
])
                   