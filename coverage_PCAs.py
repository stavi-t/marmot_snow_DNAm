#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import numpy as np
import glob
import gzip
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
#from sklearn.preprocessing import StandardScaler
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import os


# assembling dataframe for plotting coverage metrics
# set working dir
cgmap_dir = "/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/bsbolt/called/filter_10x-88x"

# initialize a list to store data from all samples
all_data = []

# glob all .CGmap files in the called directory
cgmap_files = glob.glob(os.path.join(cgmap_dir, "*.CGmap.gz"))

#loop over CGmap files and process each
for files in cgmap_files:
    # test first 5 files
    # if len(all_data)>4:
    #     break
# extract sample name from glob file name
    sample_name = os.path.basename(files).split(".")[0]

# read in and filter the CGmap file for relevant columns
    with gzip.open(files, "rt") as f:
        df = pd.read_csv(f, sep='\t', header=None, usecols=[0, 3, 7],
                     names=["Position", "Context", "Coverage"])

    # print a message indicating when each file has been read
    print(f"File {sample_name} read")

    # filter for CpG context only
    df = df[df["Context"] == "CG"]

    # add column for the sample name
    df["Sample"] = sample_name

    # Append each sample's data to all_data, which is list
    all_data.append(df)

print("All files appended.")

# make sure len=number of samples
len(all_data)

# concatenate giant dataframe (so that all_data is no longer an indexable list)
all_data = pd.concat(all_data)

###################################################################################################################
# PCA VISUALIZATIONS
###################################################################################################################
# read compiled metadata for 189 samples
df_metadata=pd.read_csv("/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/vis/metadata_189_puc19.csv")
df_metadata.head()

# read freq methylation matrix with >80% data for 165 samples - matrix then filtered to 575k, now has 10.6% missing data
# cleaned up quotation marks from file before reading in, otherwise pandas hates this
df_methyl = pd.read_csv("/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/bsbolt/filtered/575k_0.80_165_freq_cleaned.txt", delim_whitespace=True, header=0)
df_methyl.head()

# transpose and clean up matrix
df_methyl = df_methyl.T # transpose to make samples into rows, Site position into columns
df_methyl.columns = df_methyl.iloc[0] # set index row for Site
df_methyl = df_methyl.iloc[1:] # drop  the first row (Sites) because now it is the header
df_methyl.reset_index(inplace=True)
df_methyl.rename(columns={ df_methyl.columns[0]: "file_name" }, inplace = True) # add file_name column for df merge

# MERGE methylation data with metadata to continue with only 165 samples which have methylation data
df_merged = df_metadata.merge(df_methyl)
df_merged.head()
#df_merged.to_csv("/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/bsbolt/vis/methyl_metadata_165.csv")

# for plotting - drop missing CpGs from original methylation matrix (10.6%)
df_methyl_plot = df_methyl_original.dropna()
# exclude the first column containing Site
df_methyl_plot = df_methyl_plot.iloc[:,1:]

# convert df to numpy array to transpose
# once an array, no longer has column names
X = df_methyl_plot.to_numpy()
# transpose np array so samples are now rows
X = X.T

# setup PCA model with 2 PCs
qc_pca = PCA(n_components=2)

qc_pca.fit(X) # find the 2 principal components
Y = qc_pca.transform(X) # map data to principal components

# explained PC variance ratio
var_ratio = qc_pca.explained_variance_ratio_
print(f"PC1: {var_ratio[0]*100:.1f}%")
print(f"PC2: {var_ratio[1]*100:.1f}%")
#PC1: 8.1%
#PC2: 6.1%

# MANUSCRIPT PCAs
fig, axs = plt.subplots(2, 2, figsize=(8, 8))

# set metadata variables from merged 165 file
var_batch = df_merged["batch"]
var_age = df_merged["age"].astype("category")
var_year = df_merged["year"]
var_jdate = df_merged["jdate"]

# BATCH - categorical, use sbrn
sns.scatterplot(ax=axs[0, 0], x=Y[:, 0], y=Y[:, 1], hue=var_batch, palette="viridis",s=60, alpha=0.9)
axs[0, 0].set_title("a) Batch", loc="left", fontsize=14)
axs[0, 0].legend(loc="lower right", fontsize=8, frameon=True, edgecolor="grey")

# AGE - categorical, use sbrn
sns.scatterplot(ax=axs[0, 1], x=Y[:, 0], y=Y[:, 1], hue=var_age, palette="viridis", s=60, alpha=0.9)
axs[0, 1].set_title("b) Age", loc="left", fontsize=14)
axs[0, 1].legend(loc="lower right", fontsize=8, frameon=True, edgecolor="grey")

# YEAR - continuous
sc3 = axs[1, 0].scatter(Y[:, 0], Y[:, 1], c=var_year, cmap="coolwarm", s=60, alpha=1, edgecolors="white", linewidth=0.5)
axs[1, 0].set_title("c) Year", loc="left", fontsize=14)
cbar3 = plt.colorbar(sc3, ax=axs[1, 0])

# JDATE - continuous
# boolean mask filter for na vals and unreliable emergence dates past Jun 9 (from Edic et al. 2020)
mask = (~np.isnan(var_jdate)) & (var_jdate <= 160)
Y_emerg = Y[mask]
var_jdate_clean = var_jdate[mask]

sc4 = axs[1, 1].scatter(Y_emerg[:, 0], Y_emerg[:, 1],c=var_jdate_clean,cmap="coolwarm",s=60, alpha=0.9, edgecolors="white", linewidth=0.5)
axs[1, 1].set_title("d) Emergence Date", loc="left", fontsize=14)
cbar4 = plt.colorbar(sc4, ax=axs[1, 1])

# global axis labels
pc1 = qc_pca.explained_variance_ratio_[0] * 100
pc2 = qc_pca.explained_variance_ratio_[1] * 100

fig.supxlabel(f"PC1 ({pc1:.2f}%)", fontsize=14)
fig.supylabel(f"PC2 ({pc2:.2f}%)", fontsize=14)

plt.show()


df_merged["age"].values
var_jdate.values
var_jdate_clean.values

###################################################################################################################
# PCA STATS
###################################################################################################################
# PC1 and PC2 extracted directly from Y
PC1 = Y[:, 0]
PC2 = Y[:, 1]

# PC1 vs. year - continuous (with corrected vals)
pearson_coef, p_value = pearsonr(PC1, var_year)
print("Pearson corr:", pearson_coef, "p-val=", p_value)
#Pearson corr: -0.01468883290500498 p-val= 0.8514590072078582

# PC2 vs. year - continuous (with corrected vals)
pearson_coef, p_value = pearsonr(PC2, var_year)
print("Pearson corr:", pearson_coef, "p-val=", p_value)
#Pearson corr: 0.05489470573060071 p-val= 0.483738155489783
# PC2 vs. jdate - continuous (with corrected vals)
pearson_coef, p_value = pearsonr(Y_emerg[:, 1], var_jdate_clean)
print("Pearson corr:", pearson_coef, "p-val=", p_value)
#Pearson corr: 0.05041453377195458 p-val= 0.5599784441981629

# Spearman rank corr. for categorical batch
rho_batch, p_batch = spearmanr(PC1, var_batch)
print("Spearman's rank", rho_batch, "p-val", p_batch)
#Spearman's rank -0.05335906631933271 p-val 0.4960714171778301
# Spearman rank corr. for categorical batch
rho_batch, p_batch = spearmanr(PC2, df_merged["batch"])
print("Spearman's rank corr:", rho_batch, "p-val=", p_batch)
#Spearman's rank corr: 0.00843860503771015 p-val= 0.9143338601870956

# Spearman corr. for categorical age
rho_age, p_age = spearmanr(PC1, var_age)
print("Spearman's rank corr:", rho_age, "p-val=", p_age)
#Spearman's rank corr: 0.08341194697972841 p-val= 0.2868045138437879

# Spearman corr. for categorical age
rho_age, p_age = spearmanr(PC2, var_age)
print("Spearman's rank corr:", rho_age, "p-val=", p_age)
#Spearman's rank corr: 0.10466584815960285 p-val= 0.18092446751042265

