import pandas as pd 
import numpy as np
import scanpy as sc
import yaml
import scanpy as sc

#statistics
import statistical_file as sf

#math
from scipy.stats import norm

#matplotlib
from matplotlib_venn import venn2_unweighted  
from matplotlib import pyplot as plt 

#Plotly
import plotly.graph_objects as go
import plotly.express as px

#bokeh
from bokeh.plotting import output_notebook,figure, show
from bokeh.models import ColumnDataSource

#panel
import panel as pn


#scikit
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import QuantileTransformer
from sklearn.decomposition import PCA
from scipy.stats import chi2
from scipy.stats.mstats import zscore
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 
warnings.filterwarnings("ignore", category=RuntimeWarning) 

import statistical_file as sf

def get_config(file_name):
   '''This is the configuaration file'''
   with open(file_name,'r', encoding = 'UTF-8') as stream:
      config = yaml.safe_load(stream)
   return config  


config = get_config('config.yaml')
microarray_df = pd.read_table(config['microarray'])
microarray_df = microarray_df.set_index(microarray_df.iloc[:, 0])
microarray_df = microarray_df.T
microarray_df.drop(microarray_df.index[0], inplace=True)

gene_sample_data_df = pd.read_excel(config['gene_sample'], engine='odf')

lung3_meata_df = pd.read_excel(config['lung3_metadata'])
merged_lung3_data = lung3_meata_df.merge(gene_sample_data_df[['Sample_title','Sample_geo_accession']], left_on = 'title', right_on='Sample_title')
merged_lung3_data= merged_lung3_data.set_index(merged_lung3_data['Sample_geo_accession'])

final_dataset = microarray_df.merge(merged_lung3_data, left_index=True,right_index=True)
cols_to_delete = ['sample.name', 'title', 'CEL.file', 'source.location', 'organism', 'characteristics.tag.stage.mets',
       'characteristics.tag.primaryVSmets', 'characteristics.tag.grade',
       'molecule tested', 'label', 'platform', 'Sample_title']
final_dataset.drop(cols_to_delete, axis=1, inplace=True)
df_male = final_dataset[final_dataset['characteristics.tag.gender'] == 'M']
df_female = final_dataset[final_dataset['characteristics.tag.gender'] == 'F']

def statistical_test(probe_set):
   df_F = df_female[probe_set].astype('float64')
   df_M = df_female[probe_set].astype('float64')
   f_mean = df_F.mean()
   m_mean = df_M.mean()
 
   ad_test_values_F = sf.DS_AndersonDarling_test_normal(df_F)
   ad_test_values_M = sf.DS_AndersonDarling_test_normal(df_M)
   
   if (ad_test_values_F[2] and ad_test_values_M[2]) > 0.05:
      sig_test_values = sf.DS_2sample_ttest_means(df_F,df_M, alternative='two-sided')
   else:
      sig_test_values = sf.DS_2sample_MannWhitney_test_medians (df_F,df_M, alternative='two-sided')
   
   print(sig_test_values[1])
   if sig_test_values[1] < 0.05:
      sig_diff = 'Yes'
   else:
      sig_diff = 'No'
   
   stat_data = pd.DataFrame({
         'gender': ['Female','Male'],
         'mean_exp': [f_mean, m_mean],
         'ad_pvalue': [ad_test_values_F[2],ad_test_values_M[2]],
         'mann_whitney': [sig_test_values[1],sig_test_values[1]],
         'sig_diff': [sig_diff,sig_diff]
          })
   type(stat_data)
   return stat_data
   
# print(statistical_test('merck-CR749383_at'))
      

