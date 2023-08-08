# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 13:53:31 2022

@author: wadus
"""
#%%
#Libraries
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from sklearn.utils import shuffle
from sklearn.model_selection import cross_validate
from sklearn.inspection import permutation_importance
from sklearn import metrics
from matplotlib import pyplot as plt
import numpy as np
from yellowbrick.classifier import ConfusionMatrix
from yellowbrick.model_selection import FeatureImportances
#%%
# Crossvalidation
k=4
def cross_validation(model, _X, _y, _cv=k):
      '''Function to perform k Folds Cross-Validation
       Parameters
       ----------
      model: Python Class, default=None
              This is the machine learning algorithm to be used for training.
      _X: array
           This is the matrix of features.
      _y: array
           This is the target variable.
      _cv: int, default=5
          Determines the number of folds for cross-validation.
       Returns
       -------
       The function returns a dictionary containing the metrics 'accuracy', 'precision',
       'recall', 'f1' for both training set and vlidation set.
      '''
      _scoring = ['accuracy', 'precision', 'recall', 'f1']
      results = cross_validate(estimator=model,
                               X=_X,
                               y=_y,
                               cv=_cv,
                               scoring=_scoring,
                               return_train_score=True)
      
      return {
              "Mean Validation Accuracy": results['test_accuracy'].mean()*100,
              "Mean Validation Precision": results['test_precision'].mean(),
              "Mean Validation Recall": results['test_recall'].mean(),
              "Mean Validation F1 Score": results['test_f1'].mean()
              }
#%%
#Load df
mt="C:\\Users\\wadus\\OneDrive\\Desktop\\Studium\\Masterarbeit\\Projects\\Final_Data\\ML_model\\ml_hvc.csv"

df = pd.read_csv(mt)
df = shuffle(df)
X=df.iloc[:,2:-1]
y=df.iloc[:,-1]
# split X and y into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=42,stratify=y)
# instantiate the model (using the default parameters)
LOG_model = LogisticRegression()
KNN_model = KNeighborsClassifier(n_neighbors=5)
GNB_model = GaussianNB()
DTC_model = DecisionTreeClassifier()
GBC_model = GradientBoostingClassifier()
# fit the model with data
LOG_model.fit(X_train, y_train)
KNN_model.fit(X_train, y_train)
GNB_model.fit(X_train,y_train)
DTC_model.fit(X_train,y_train)
GBC_model.fit(X_train,y_train)

#Calculation of result
LOG_result = cross_validation(LOG_model, X, y,k)
KNN_result = cross_validation(KNN_model, X, y,k)
GNB_result = cross_validation(GNB_model, X, y,k)
DTC_result = cross_validation(DTC_model, X, y,k)
GBC_result = cross_validation(GBC_model, X, y,k)
#%%
#Printing
print("LOG results \n",LOG_result)
print("KNN results \n",KNN_result)
print("GNB results \n",GNB_result)
print("DTC results \n",DTC_result)
print("GBC results \n",GBC_result)
#%%
#set up plotting area
plt.figure(0).clf()

#fit logistic regression model and plot ROC curve
model = LogisticRegression()
model.fit(X_train, y_train)
y_pred = model.predict_proba(X_test)[:, 1]
fpr, tpr, _ = metrics.roc_curve(y_test, y_pred)
auc = round(metrics.roc_auc_score(y_test, y_pred), 4)
plt.plot(fpr,tpr,label="LOG, AUC="+str(auc))

#fit knn model and plot ROC curve
model = KNeighborsClassifier()
model.fit(X_train, y_train)
y_pred = model.predict_proba(X_test)[:, 1]
fpr, tpr, _ = metrics.roc_curve(y_test, y_pred)
auc = round(metrics.roc_auc_score(y_test, y_pred), 4)
plt.plot(fpr,tpr,label="KNN, AUC="+str(auc))

#fit gnb model and plot ROC curve
model = GaussianNB()
model.fit(X_train, y_train)
y_pred = model.predict_proba(X_test)[:, 1]
fpr, tpr, _ = metrics.roc_curve(y_test, y_pred)
auc = round(metrics.roc_auc_score(y_test, y_pred), 4)
plt.plot(fpr,tpr,label="GNB, AUC="+str(auc))

#fit dtc model and plot ROC curve
model = DecisionTreeClassifier()
model.fit(X_train, y_train)
y_pred = model.predict_proba(X_test)[:, 1]
fpr, tpr, _ = metrics.roc_curve(y_test, y_pred)
auc = round(metrics.roc_auc_score(y_test, y_pred), 4)
plt.plot(fpr,tpr,label="DTC, AUC="+str(auc))

#fit gradient boosted model and plot ROC curve
model = GradientBoostingClassifier()
model.fit(X_train, y_train)
y_pred = model.predict_proba(X_test)[:, 1]
fpr, tpr, _ = metrics.roc_curve(y_test, y_pred)
auc = round(metrics.roc_auc_score(y_test, y_pred), 4)
plt.plot(fpr,tpr,label="GBC, AUC="+str(auc))

#add legend
plt.legend()
plt.savefig("C:\\Users\\wadus\\OneDrive\\Desktop\\Studium\\Masterarbeit\\Projects\\Final_Plots\\ML model\\hvc_roc.svg")
#%%
model = GradientBoostingClassifier()
viz = FeatureImportances(model)
#del X_train['p201_250']
#del X_train['p181_200']
#del X_train['p151_180']
viz.fit(X_train, y_train)
viz.show("C:\\Users\\wadus\\OneDrive\\Desktop\\Studium\\Masterarbeit\\Projects\\Final_Plots\\ML model\\feature_importance.svg")
#%%
# Crossvalidation
k=4
def multi_cross_validation(model, _X, _y, _cv=k):
      '''Function to perform k Folds Cross-Validation
       Parameters
       ----------
      model: Python Class, default=None
              This is the machine learning algorithm to be used for training.
      _X: array
           This is the matrix of features.
      _y: array
           This is the target variable.
      _cv: int, default=5
          Determines the number of folds for cross-validation.
       Returns
       -------
       The function returns a dictionary containing the metrics 'accuracy', 'precision',
       'recall', 'f1' for both training set and vlidation set.
      '''
      _scoring = ['accuracy', 'precision_macro', 'recall_macro', 'f1_macro']
      results = cross_validate(estimator=model,
                               X=_X,
                               y=_y,
                               cv=_cv,
                               scoring=_scoring,
                               return_train_score=True)
      
      return {
              "Mean Validation Accuracy": results['test_accuracy'].mean()*100,
              "Mean Validation Precision": results['test_precision_macro'].mean(),
              "Mean Validation Recall": results['test_recall_macro'].mean(),
              "Mean Validation F1 Score": results['test_f1_macro'].mean()
              }
#%%
#Load df
mt="C:\\Users\\wadus\\OneDrive\\Desktop\\Studium\\Masterarbeit\\Projects\\Final_Data\\ML_model\\ml_too_2.csv"

df = pd.read_csv(mt)
df = shuffle(df)
X=df.iloc[:,2:-1]
y=df.iloc[:,-1]
# split X and y into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=42,stratify=y)
# instantiate the model (using the default parameters)
LOG_model = LogisticRegression()
KNN_model = KNeighborsClassifier(n_neighbors=5)
GNB_model = GaussianNB()
DTC_model = DecisionTreeClassifier()
GBC_model = GradientBoostingClassifier()
# fit the model with data
LOG_model.fit(X_train, y_train)
KNN_model.fit(X_train, y_train)
GNB_model.fit(X_train,y_train)
DTC_model.fit(X_train,y_train)
GBC_model.fit(X_train,y_train)
#%%
#Calculation of result
LOG_result = multi_cross_validation(LOG_model, X, y,k)
KNN_result = multi_cross_validation(KNN_model, X, y,k)
GNB_result = multi_cross_validation(GNB_model, X, y,k)
DTC_result = multi_cross_validation(DTC_model, X, y,k)
GBC_result = multi_cross_validation(GBC_model, X, y,k)
#%%
#Printing
print("LOG results \n",LOG_result)
print("KNN results \n",KNN_result)
print("GNB results \n",GNB_result)
print("DTC results \n",DTC_result)
print("GBC results \n",GBC_result)
#%%
model = GradientBoostingClassifier()
viz = FeatureImportances(model)
viz.fit(X_train, y_train)
viz.show("C:\\Users\\wadus\\OneDrive\\Desktop\\Studium\\Masterarbeit\\Projects\\Final_Plots\\ML model\\feature_importance_too.svg")
#%%
# ROCAUC lOG
visualizer = ROCAUC(LOG_model, classes=["Healthy", "Lung", "Breast","Colorectal","Pancreas","Bile duct","Stomach","Ovarian"])

visualizer.fit(X_train, y_train)
visualizer.score(X_test, y_test)        
visualizer.show(outpath="C:\\Users\\wadus\\OneDrive\\Desktop\\Studium\\Masterarbeit\\Projects\\Final_Plots\\ML model\\log_too_roc.svg")
#%%
# The ConfusionMatrix visualizer taxes a model
cm = ConfusionMatrix(LOG_model, classes=["Healthy", "Lung", "Breast","Colorectal","Pancreas","Bile duct","Stomach","Ovarian"])

# Fit fits the passed model. This is unnecessary if you pass the visualizer a pre-fitted model
cm.fit(X_train, y_train)

# To create the ConfusionMatrix, we need some test data. Score runs predict() on the data
# and then creates the confusion_matrix from scikit-learn.
cm.score(X_test, y_test)

# How did we do?
cm.show(outpath="C:\\Users\\wadus\\OneDrive\\Desktop\\Studium\\Masterarbeit\\Projects\\Final_Plots\\ML model\\log_too_cm.svg")
#%%                   
# ROCAUC KNN
visualizer = ROCAUC(KNN_model, classes=["Healthy", "Lung", "Breast","Colorectal","Pancreas","Bile duct","Stomach","Ovarian"])

visualizer.fit(X_train, y_train)
visualizer.score(X_test, y_test)        
visualizer.show(outpath="C:\\Users\\wadus\\OneDrive\\Desktop\\Studium\\Masterarbeit\\Projects\\Final_Plots\\ML model\\knn_too_roc.svg") 
#%%
# The ConfusionMatrix visualizer taxes a model
cm = ConfusionMatrix(KNN_model, classes=["Healthy", "Lung", "Breast","Colorectal","Pancreas","Bile duct","Stomach","Ovarian"])

# Fit fits the passed model. This is unnecessary if you pass the visualizer a pre-fitted model
cm.fit(X_train, y_train)

# To create the ConfusionMatrix, we need some test data. Score runs predict() on the data
# and then creates the confusion_matrix from scikit-learn.
cm.score(X_test, y_test)

# How did we do?
cm.show(outpath="C:\\Users\\wadus\\OneDrive\\Desktop\\Studium\\Masterarbeit\\Projects\\Final_Plots\\ML model\\knn_too_cm.svg")
#%%
# ROCAUC GNB
visualizer = ROCAUC(GNB_model, classes=["Healthy", "Lung", "Breast","Colorectal","Pancreas","Bile duct","Stomach","Ovarian"])

visualizer.fit(X_train, y_train)
visualizer.score(X_test, y_test)        
visualizer.show(outpath="C:\\Users\\wadus\\OneDrive\\Desktop\\Studium\\Masterarbeit\\Projects\\Final_Plots\\ML model\\gnb_too_roc.svg")
#%%
# The ConfusionMatrix visualizer taxes a model
cm = ConfusionMatrix(GNB_model, classes=["Healthy", "Lung", "Breast","Colorectal","Pancreas","Bile duct","Stomach","Ovarian"])

# Fit fits the passed model. This is unnecessary if you pass the visualizer a pre-fitted model
cm.fit(X_train, y_train)

# To create the ConfusionMatrix, we need some test data. Score runs predict() on the data
# and then creates the confusion_matrix from scikit-learn.
cm.score(X_test, y_test)

# How did we do?
cm.show(outpath="C:\\Users\\wadus\\OneDrive\\Desktop\\Studium\\Masterarbeit\\Projects\\Final_Plots\\ML model\\gnb_too_cm.svg") 
#%%
# ROCAUC DTC
visualizer = ROCAUC(DTC_model, classes=["Healthy", "Lung", "Breast","Colorectal","Pancreas","Bile duct","Stomach","Ovarian"])

visualizer.fit(X_train, y_train)
visualizer.score(X_test, y_test)        
visualizer.show(outpath="C:\\Users\\wadus\\OneDrive\\Desktop\\Studium\\Masterarbeit\\Projects\\Final_Plots\\ML model\\dtc_too_roc.svg")
#%%
# The ConfusionMatrix visualizer taxes a model
cm = ConfusionMatrix(DTC_model, classes=["Healthy", "Lung", "Breast","Colorectal","Pancreas","Bile duct","Stomach","Ovarian"])

# Fit fits the passed model. This is unnecessary if you pass the visualizer a pre-fitted model
cm.fit(X_train, y_train)

# To create the ConfusionMatrix, we need some test data. Score runs predict() on the data
# and then creates the confusion_matrix from scikit-learn.
cm.score(X_test, y_test)

# How did we do?
cm.show(outpath="C:\\Users\\wadus\\OneDrive\\Desktop\\Studium\\Masterarbeit\\Projects\\Final_Plots\\ML model\\dtc_too_cm.svg") 
#%%
# ROCAUC GBC
visualizer = ROCAUC(GBC_model, classes=["Healthy", "Lung", "Breast","Colorectal","Pancreas","Bile duct","Stomach","Ovarian"])

visualizer.fit(X_train, y_train)
visualizer.score(X_test, y_test)        
visualizer.show(outpath="C:\\Users\\wadus\\OneDrive\\Desktop\\Studium\\Masterarbeit\\Projects\\Final_Plots\\ML model\\gbc_too_roc.svg") 
#%%
# The ConfusionMatrix visualizer taxes a model
cm = ConfusionMatrix(GBC_model, classes=["Healthy", "Lung", "Breast","Colorectal","Pancreas","Bile duct","Stomach","Ovarian"])

# Fit fits the passed model. This is unnecessary if you pass the visualizer a pre-fitted model
cm.fit(X_train, y_train)

# To create the ConfusionMatrix, we need some test data. Score runs predict() on the data
# and then creates the confusion_matrix from scikit-learn.
cm.score(X_test, y_test)

# How did we do?
cm.show(outpath="C:\\Users\\wadus\\OneDrive\\Desktop\\Studium\\Masterarbeit\\Projects\\Final_Plots\\ML model\\gbc_too_cm.svg")