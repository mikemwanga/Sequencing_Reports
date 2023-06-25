
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
import dash_ag_grid as dag
load_figure_template("pulse")

PATH = pathlib.Path(__file__).parent
DATA_PATH = PATH.joinpath("./data").resolve()
# #cerulean,flatly,journal,litera,pulse,sandstone,minty


col_title = "text-center text-black fw-normal"
style_text ={"font-size":12,"text-align":"center"}

col_title = "text-center text-black fw-normal"
markercolor = "#8B0000"
margin = dict(l=10, r=10, t=5, b=5) #'#543005'
color_patterns = ['#1b9e77','#d95f02','#7570b3','magenta','#00ffef','#40004b','#ff0000','#053061','#cae00d','#4169e1',"#67001f"]
pcolor = "rgba(0,0,0,0)"
pcolor_white = "white"
gridcolor="lightgray"
pcolor_home = "#E6E6E6"
cardbody_style = {"background-color":pcolor}
style_text ={"font-size":14,"text-align":"center"}