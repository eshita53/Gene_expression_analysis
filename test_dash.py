import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import dash
from dash import dcc
from dash.dependencies import Input, Output
import load_df as data 
import plotly.graph_objects as go
from dash import Dash, dcc, html, Input, Output,dash_table,ctx
import dash_bootstrap_components as dbc
from dash_bootstrap_templates import load_figure_template


# app = Dash(__name__)
dbc_css = ("https://cdn.jsdelivr.net/gh/AnnMarieW/dash-bootstrap-templates@V1.0.2/dbc.min.css")
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP, dbc_css])

DATA_TABLE_COLUMNS = [
    {
        "id": "gender",
        "name": "Gender"
    },
     {
        "id": "mean_exp",
        "name": "Mean",
        "type": "numeric"
    },
    {
        "id": "ad_pvalue",
        "name": "Anderson Test (p_value)",
        "type": "numeric"
    },
    { "id": "mann_whitney", 
     "name": "MannWhitney Test (p_value)"},
    
    { "id": "sig_diff", 
     "name": "Significant Diff"}
]

DATA_TABLE_STYLE = {
    "style_header": {
        "color": "white",
        "backgroundColor": "#799DBF",
        "fontWeight": "bold",
    },
}
def format_table_row(row):

    row['mean_exp'] = round(row['mean_exp'], 4)
    row['ad_pvalue'] = round(row['ad_pvalue'], 4)
    
    return row

# app.layout = html.Div( [
app.layout = dbc.Container(
    [
    html.H1("Single Probset Analysis on Male and Female Patients", style={'text-align': 'center'}),

    dcc.Dropdown(
        id = 'probe_set',
        options = [
            {'label': probe_set, 'value': probe_set}
            for probe_set in data.final_dataset.columns[:50]
        ], 
        clearable = False,
        searchable = False,
        className = 'dropdown', style={'fontSize': "24px",'textAlign': 'center'},
    ),
    html.Div(id='output_container', children=[]),
    html.Br(),
    html.Div(
        [
            dbc.Row(
                
                [
                    dbc.Col( dcc.Graph(id='box_plot', figure={}, className="border"), lg=6),
                    dbc.Col( dcc.Graph(id='histogram', figure={}, className="border"), lg=6),
                ],
                className="mt-4",
            ),

            dbc.Row(
                [
                    dbc.Col( dcc.Graph(id='violinplot', figure={}, className="border"), lg=6),
                    dbc.Col( 
                    dash_table.DataTable(
                    id="user_datatable",
                    sort_action="native",
                    style_header=DATA_TABLE_STYLE.get("style_header"),
                    ),lg=6),

                ],
                className="mt-5",
            )

        ]
    )    

    ],
    fluid=True,
)


# Connect the Plotly graphs with Dash Components
@app.callback(
    [
     Output(component_id='box_plot', component_property='figure'),
     Output(component_id='histogram', component_property='figure'),
     Output(component_id='user_datatable', component_property='data'),
     Output(component_id='user_datatable', component_property='columns'),
     Output(component_id='violinplot', component_property='figure'), 
   ],
    [Input(component_id='probe_set', component_property='value')]
)
def update_graph(probe_set):

    boxplot = go.Figure(
        data =[
            (go.Box(y=data.df_female[probe_set], name = 'Female', hovertemplate=None)),
            (go.Box(y=data.df_male[probe_set], name = 'Male', hovertemplate=None)),    
        ]
        
    )
    boxplot.update_layout( boxmode='group', 
                boxgroupgap=0.2,
                title='Probe Set profiles of Female and Male Patients', 
                yaxis=dict(title="{} values ".format(probe_set)),
                hovermode="x", 
                font = dict(size = 14, family="Arial Black"),
                margin=dict(l=100, r=100, b=100, t=80))
    #histogram
    fig2 = go.Figure(
        data=[
        go.Histogram(
            x=data.df_female[probe_set], name ='Female'
            ),
        go.Histogram(
            x=data.df_female[probe_set], name ='Male'
            )
        ])
    fig2.update_layout(barmode='stack',title='Female and Male Patients probset expression' )
    fig2.update_traces(opacity=0.75)
    fig2.update_xaxes(title_text='Value of expression profiles')
    fig2.update_yaxes(title_text='Frequency')
    
    #violin
    fig = go.Figure()
    fig.add_trace(go.Violin(y=data.df_female[probe_set], box_visible=True, line_color='black',
                                meanline_visible=True, fillcolor='lightseagreen', opacity=0.6,
                                x0='Female', name = 'Female'))
    fig.add_trace(go.Violin(y= data.df_female[probe_set], box_visible=True, line_color='black',
                                meanline_visible=True, fillcolor='dodgerblue', opacity=0.6,
                                x0='Male', name = 'Male'))
    fig.update_layout(title='Violin plot of Female and Male Patients gene expression profiles', 
                    yaxis=dict(title='Gene Expression Profiles'),
                    hovermode="x", 
                    font = dict(size = 14, family="Arial Black"),
                    margin=dict(l=100, r=100, b=100, t=80),
                  legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1))


    
    
    stat_data_df = data.statistical_test(probe_set)
    updated_table = stat_data_df.to_dict('records')
    
    updated_table = [format_table_row(row) for row in updated_table]



    return  boxplot, fig2, updated_table, DATA_TABLE_COLUMNS,fig


if __name__ == '__main__':
    app.run_server()