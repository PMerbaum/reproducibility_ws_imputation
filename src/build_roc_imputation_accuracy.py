#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 14:37:18 2024

@author: pmerbaum

This code is an example for short tandem repeat (STR) imputation accuracy analysis.
Allele lengths may be cruisial in some cases. For patients with amyotrophic lateral sclerosis,
a C9ortf72 gene length > 30 repeats means shorter survival and more aggressive disease.
We want to establish a model that will predict (impute) these lengths in the most accurate way. 
To test our model, we took samples with known STR lengths, masked C9ortf72 gene region, 
reimputed it and now we need to check, how accurate our imputation is.

We have columns in df 'data':
    
data['imputed'] - calculated length on 2 alleles (format: l1/l2);
data['true'] - known length on 2 alleles (format: l1/l2);
data['ext_dosage'] - probability that this allele is extended (>30 STR) calculated by the software (formar: float);
data['GWAS_SNP_dosage'] - probability that the most significant SNP in our GWAS is present (formar: float). We won't use it for now.
"""
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
from sklearn.metrics import ConfusionMatrixDisplay
from sklearn.metrics import roc_curve, auc

data = pd.read_csv("data/raw/fake_data_for_ws.txt", sep="\t")

# swap allele lengths to short/long independed from phasing
def swap_values(data, column_name):
    split_values = data[column_name].str.split("/")
    swapped_values = split_values.apply(
        lambda x: f"{x[1]}/{x[0]}" if int(x[0]) > int(x[1]) else "/".join(x)
    )
    data.loc[:, column_name] = swapped_values

swap_values(data, 'imputed')
swap_values(data, 'true')

# creating columns for short and long lengths
def split_short_long(data, columns):
    for column in columns:
        split_values = data[column].str.split("/")
        split_values = split_values.apply(lambda x: [int(e) for e in x])

        data.loc[:, f"{column}_short"] = split_values.str[0].astype(int)
        data.loc[:, f"{column}_long"] = split_values.str[1].astype(int)

split_short_long(data, ['imputed', 'true'])

# now, to create a confusion matrix (2x2) for correct and incorrect imputed lengths,
# we need to transform longest lengths into binary values.

imp = data['imputed_long'].values
true = data['true_long'].values
imp_lst = imp.tolist()
true_lst = true.tolist()

threshold = 30
binary_values_imp = [1 if x >= threshold else 0 for x in imp_lst]
binary_values_true = [1 if x >= threshold else 0 for x in true_lst]

cm = confusion_matrix(binary_values_true, binary_values_imp)
disp = ConfusionMatrixDisplay(confusion_matrix=cm)
disp.plot()
plt.savefig("confusion_matrix_imputation.png")

# build a ROC curve for different dosages

true = data.iloc[:, 3].values
true_m = true.reshape(-1, 1)
binary_values_eh = [1 if x >= threshold else 0 for x in true_m]
binary_values_imp = data.iloc[:, 6].values

fpr, tpr, thresholds = roc_curve(binary_values_true, binary_values_imp)
roc_auc = auc(fpr, tpr)

plt.clf()
plt.figure()
plt.plot(
    fpr,
    tpr,
    color="green",
    lw=2,
    label="ROC curve for imputation (AUC = %0.2f)" % roc_auc,
)


plt.plot([0, 1], [0, 1], color="navy", lw=2, linestyle="--")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("false positive rate")
plt.ylabel("true positive rate ")
plt.title("ROC for STR imputation - dosage")
plt.legend(loc="lower right")

plt.show()
plt.savefig("roc_imputation_accuracy.png")
