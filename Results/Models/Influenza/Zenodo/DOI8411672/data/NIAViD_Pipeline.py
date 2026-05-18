# -*- coding: utf-8 -*-
"""
@author: alphaforna
"""

# load required packages
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.image as mpimg
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.ensemble import IsolationForest
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import mean_squared_error
from sklearn import model_selection
from sklearn import metrics
import eli5
from eli5.sklearn import PermutationImportance

import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerBase
from matplotlib.lines import Line2D
from sklearn.svm import SVC
from sklearn.svm import OneClassSVM
from sklearn.neighbors import LocalOutlierFactor
from sklearn.inspection import permutation_importance
from sklearn.metrics import f1_score
from sklearn.metrics import recall_score
from sklearn.metrics import precision_score
from sklearn.metrics import roc_auc_score

from sklearn.utils import shuffle
from matplotlib.legend_handler import HandlerBase
from matplotlib.lines import Line2D
from math import sqrt
from pyod.models.ecod import ECOD
from statsmodels.distributions.empirical_distribution import ECDF
from statsmodels.stats.proportion import proportion_confint

import hdbscan

import seaborn as sns
import sklearn.datasets as data

import warnings

# set outlier detection algorithm
# outlier_method = "SVM"
# outlier_method = "iForest"
# outlier_method = "ECOD"
outlier_method = "HDBSCAN"
def setLabel():
  if outlier_method == "iForest" or outlier_method == "SVM":
    label = [-1,1]
  if outlier_method == "ECOD" or outlier_method == "HDBSCAN":
    label = [1,0]
  return label
transLabel = setLabel()

# outlier threshold for HDBSCAN
threshold = 0.5

# all features in dataset
FEATURES=['distRoot', 'branch_length', 'hydrophobicity', 'charge', 'boman', 'instability', 'isoelectric_point'] 
# Load data and data wrangling
# Load All the sites
file_path="/path/to/sequence_physioproperties.csv"
physio_data=pd.read_csv(file_path)

# Get numerical features and scale them to zero mean and unit variance
ss=StandardScaler()
to_scale=physio_data[FEATURES]
scaled_df=pd.DataFrame(ss.fit_transform(to_scale), 
                       columns=to_scale.columns)

scaled_df=physio_data[FEATURES]

# scale to [0,1] 
for i in range(2,7):
    scaled_df.iloc[:,i] = (scaled_df.iloc[:,i]-np.min(scaled_df.iloc[:,i]))/(np.max(scaled_df.iloc[:,i])-np.min(scaled_df.iloc[:,i]))
    

# Add the categorical variable back in
scaled_df['cluster'] = physio_data['cluster']

# Add the sequence IDs back too.
scaled_df['seq_id'] = physio_data['seq_id']
##################################################################
year=[x.split('/')[2] for x in physio_data['seq_id'].to_list()]
for i in range(0, len(year)):
    if int(year[i]) > 67:
        year[i]=int('19'+year[i])
    else:
        year[i]=int('20'+year[i])

physio_data = scaled_df
physio_data['year']=year

def cluster_data(dat):
    """(0: inliers, 1: outliers) for ECOD
       (1: inliers, -1: outliers) for ECOD
    """
    dat.reset_index(inplace=True, drop =True)
    uniqueValues, indicesList = np.unique(dat.cluster, return_index=True)
    dat['cluster2'] = np.repeat(transLabel[1],len(dat.cluster)).copy() # non-transition label 1
    dat.loc[ :,'cluster2'][indicesList] =transLabel[0]    # transition label-1
    dat['cluster2'] = dat['cluster2'].astype('category')

def cluster_data_null(dat): 
    dty = [np.random.choice(transLabel,p=[round(11/dat.shape[0],5),1-round(11/dat.shape[0],5)]) for toss in range(dat.shape[0])]
    dat1 =dat.assign(cluster2=dty)
    dat1['cluster2'] = dat1['cluster2'].astype('category')
    return dat1

### Helper functions for misc use
# Get unique items in a list
def unique(list):
    unique_items=[]
    for i in list:
        if i not in unique_items:
            unique_items.append(i)
        else:
            pass
    return unique_items
# Constants used throughout cells

# Colors for plotting by cluster
CLUSTER_COLORS={'HK68':'dodgerblue',
                'EN72':'magenta',
                'VI75':'olive', 
                'TX77':'aqua', 
                'BK79':'red', 
                'SI87':'gold', 
                'BE89':'grey',
                'BE92':'salmon', 
                'WU95':'teal', 
                'SY97':'maroon', 
                'FU02':'palegreen'}

# Features excluding phylogenetic measurements
NO_PHYLO_FEATURES=['hydrophobicity','charge','boman','instability','isoelectric_point']
# Function to calculate mathew's correlation coefficient 
def mcc(tp, fp, tn, fn):
    # https://stackoverflow.com/a/56875660/992687
    x = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    return ((tp * tn) - (fp * fn)) / sqrt(x)
  
# Sampling Strategy
np.random.seed(15)
# Create the training (80%) an testing (20%) datasets
df_train = physio_data.groupby('cluster',as_index = False,group_keys=False).apply(lambda s: s.sample(frac=0.8))
# Validation data- Not used in buidling the model and testing it
df_test = physio_data[~physio_data.index.isin(df_train.index)]
df_train.sort_values(by=['year'], inplace=True)
df_test.sort_values(by=['year'], inplace=True)
# Transition variable reflecting true antigenic transition
cluster_data(df_train)
cluster_data(df_test)

# Training Pipeline 
''' Training Phase of the pipeline'''

train_data =df_train[df_train['year'] > 1972]
def trainFunct(dataset):
  np.random.seed(15)
  initial_data = df_train[df_train['year'] <= 1972]
  current_data = initial_data.copy()
  X = dataset[NO_PHYLO_FEATURES]
  total= []
  for y in range(0,len(dataset)):
    row = pd.DataFrame(X.iloc[y]).T.values
    next_data =dataset.iloc[y]
    X_train= current_data[NO_PHYLO_FEATURES].values
    if outlier_method == "SVM":
        iSVM  = OneClassSVM(gamma='scale')
        clf = iSVM.fit(X_train)
        y_pred_train=clf.predict(row) 
        y_prob_train=clf.decision_function(row)
    if outlier_method == "iForest":
      iForest= IsolationForest(n_estimators =100, max_samples ='auto', contamination = 'auto', random_state=0)
      clf = iForest.fit(X_train)
      y_pred_train=clf.predict(row) 
      y_prob_train=clf.decision_function(row)
    if outlier_method == "ECOD":
      clf=ECOD()
      clf.fit(X_train)
      y_pred_train=clf.predict(row) 
      y_prob_train=clf.decision_function(row)
    if outlier_method == "HDBSCAN":
        clf = hdbscan.HDBSCAN(min_cluster_size=2, gen_min_span_tree=True)
        clf.fit(X_train)
        if clf.outlier_scores_[len(current_data) - 1] > threshold:
            y_pred_train = 1
        else:
            y_pred_train = 0
        y_prob_train = clf.outlier_scores_[len(current_data) - 1] 
        if np.isnan(y_prob_train):
            y_prob_train = 1
    total.append((y_pred_train,y_prob_train))
    current_data = current_data.append(next_data)
  df1 =pd.DataFrame(total,columns=('Prediction','DecisionScore')) 
  return df1,clf

def perfMeasures(dataset,df1,clf):
  #Create categorical variables from the predictions and observations
  dataset['cluster2'] = dataset['cluster2'].astype('category')
  y_pred=np.array(df1.Prediction).astype(int)
  dataset['y_pred'] =y_pred
  dataset['y_pred'] = dataset['y_pred'].astype('category')
  y_prob=np.array(df1.DecisionScore).astype(np.float32)
  #Confusion matrix
  conf_matrix_train = confusion_matrix(dataset.cluster2,dataset.y_pred,labels=transLabel)
  # Calculate the Performance measures
  TN = conf_matrix_train[1][1]
  TP = conf_matrix_train[0][0]
  FN = conf_matrix_train[0][1]
  FP = conf_matrix_train[1][0]
  #calculate accuracy
  conf_accuracy = (float (TP+TN) / float(TP + TN + FP + FN))
  #calculate mis-classification
  conf_misclassification = 1-conf_accuracy
  
  #calculate the sensitivity
  conf_sensitivity = (TP / float(TP + FN))
  #calculate the specificity
  conf_specificity = (TN / float(TN + FP))
  # calculate recall
  conf_recall = (TP / float(TP + FN))
  # Precision
  conf_precision = (TP / float(TP + FP))
  # FI score
  conf_f1 = 2 * ((conf_precision * conf_sensitivity) / (conf_precision + conf_sensitivity))
  # AUC
  roc_auc = metrics.roc_auc_score(dataset.cluster2,y_prob)
  # MCC
  conf_mcc = mcc(tp=TP, fp=FP, tn=TN, fn=FN)
  
  #TPR vs FPR
  fpr, tpr, _ = metrics.roc_curve(dataset.cluster2, y_prob)
  # POS_LABEL??????
  # Precision recall
  prec, recall, _ = metrics.precision_recall_curve(dataset.cluster2, y_prob,pos_label=transLabel[0])
  
  #Training performance
  print(f'Accuracy: {round(conf_accuracy,5)}')
  print(f'Mis-Classification: {round(conf_misclassification,5)}')
  print(f'Sensitivity: {round(conf_sensitivity,5)}')
  print(f'Specificity: {round(conf_specificity,5)}')
  print(f'Precision: {round(conf_precision,5)}')
  print(f'f_1 Score: {round(conf_f1,5)}')
  print(f'ROC_AUC: {round(roc_auc,5)}')
  
  perfArray = [conf_accuracy,conf_misclassification,conf_sensitivity,conf_recall,conf_precision,conf_f1,roc_auc,conf_mcc,fpr,tpr,prec,recall,conf_matrix_train]
  return perfArray

def testFunct(dataset):
  np.random.seed(15)
  initial_data = df_train
  current_data = initial_data.copy()
  X = dataset[NO_PHYLO_FEATURES]
  total= []
  for y in range(0,len(dataset)):
    if dataset.iloc[y].cluster2==transLabel[1]:#Non-transition
       current_data = current_data.copy()
       current_data.sort_values(by=['year'], inplace=True)
    else:
       train_clusters=unique(current_data['cluster'].to_list()) # Get training clusters
       test_cluster=dataset.iloc[y].cluster # Get testing clusters
       train_clusters.remove(test_cluster) # Remove testing cluster from training clusters
       current_data = current_data.loc[current_data['cluster'].isin(train_clusters)] #Subset data with testing cluster removed
       current_data.sort_values(by=['year'], inplace=True)
    row = pd.DataFrame(X.iloc[y]).T.values
    next_data =dataset.iloc[y]
    ########################################################
    X_test = current_data[NO_PHYLO_FEATURES].values
    if outlier_method == "SVM":
        iSVM  = OneClassSVM(gamma='scale')
        clf = iSVM.fit(X_test)
        y_pred_test=clf.predict(row) 
        y_prob_test=clf.decision_function(row)
    if outlier_method == "iForest":
      iForest= IsolationForest(n_estimators =100, max_samples ='auto', contamination = 'auto', random_state=0)
      clf =iForest.fit(X_test)
      y_pred_test=clf.predict(row)
      y_prob_test=clf.decision_function(row)
    if outlier_method == "ECOD":
      clf =ECOD()
      clf.fit(X_test)
      y_pred_test=clf.predict(row)
      y_prob_test=clf.decision_function(row)
    if outlier_method == "HDBSCAN":
        clf = hdbscan.HDBSCAN(min_cluster_size=2, gen_min_span_tree=True)
        clf.fit(X_test)
        if clf.outlier_scores_[len(current_data) - 1] > threshold:
            y_pred_test = 1
        else:
            y_pred_test = 0
        y_prob_test = clf.outlier_scores_[len(current_data) - 1] 
        if np.isnan(y_prob_test):
            y_prob_test = 1
    ##################################################################

    total.append((y_pred_test,y_prob_test))

    current_data = current_data.append(next_data)

  df1 =pd.DataFrame(total,columns=('Prediction','DecisionScore'))
  return df1,clf

train_results = trainFunct(train_data)
df_train_result = train_results[0]
clf_train = train_results[1]
perfTrain = perfMeasures(train_data,df_train_result,clf_train)

# # with null data
# df_train_null  = cluster_data_null(df_train)
# df_test_null  = cluster_data_null(df_test)
# train_data =df_train_null[df_train_null['year'] > 1972]
# null_results = trainFunct(train_data)
# df_null_if = null_results[0]
# clf_null_if = null_results[1]
# perfNull = perfMeasures(train_data,df_null_if,clf_null_if)

df_test_sh = df_test.sample(frac = 1)
test_data =df_test_sh.copy()
test_results = testFunct(test_data)
df_test_result = test_results[0]
clf_test = test_results[1]
perfTest = perfMeasures(test_data,df_test_result,clf_test)

# plots
plt.plot(perfTrain[8],perfTrain[9])
plt.show()

plt.plot(perfTrain[11],perfTrain[10])
plt.show()

plt.plot(perfTest[8],perfTest[9])
plt.show()

plt.plot(perfTest[11],perfTest[10])
plt.show()
