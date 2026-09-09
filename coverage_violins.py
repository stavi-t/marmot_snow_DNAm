import pandas as pd
import seaborn as sns
import numpy as np
import glob
import matplotlib.pyplot as plt
import os

cd /home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/bsbolt/called/filter_10x-88x

# assembling dataframe for plotting coverage metrics - for all files in the called cov directory
# set working dirr
cgmap_dir = "/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/bsbolt/called/filter_10x-88x"

# initialize a list to store data from all samples
all_data = []

# glob all .CGmap files in the called directory
cgmap_files = glob.glob(os.path.join(cgmap_dir, "*.CGmap"))

# loop over CGmap files and process each
for files in cgmap_files:
    # test first 5 files
    # if len(all_data)>4:
    #     break
    # extract sample name from glob file name
    sample_name = os.path.basename(files).split(".")[0]

    # read in and filter the CGmap file for relevant columns
    df = pd.read_csv(files, sep='\t', header=None, usecols=[0, 3, 7],
                     names=["Position", "Context", "Coverage"])
    # Print a message indicating the file has been read
    print(f"File {sample_name} read")

    # filter for CpG context only
    df = df[df["Context"] == "CG"]

    # add column for the sample name
    df["Sample"] = sample_name

    # Append each sample's data to all_data, which is list
    all_data.append(df)

print("All files appended.")

# assembling dataframe for plotting coverage metrics - just for retained n=165 samples
# set working dirr
cgmap_dir = "/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/bsbolt/called/filter_10x-88x"
file_list_path = "/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/bsbolt/CGmaps_omit.txt"

# read full paths from file
with open(file_list_path, "r") as f:
    cgmap_files = [line.strip() for line in f if line.strip()]

# initialize a list to store data from all samples
subset_data = []

# loop over CGmap files and process each
for file_path in cgmap_files:
    sample_name = os.path.basename(file_path).split(".")[0]

    # read ZIPPED CGmap files with only the relevant columns
    df = pd.read_csv(file_path, sep='\t', header=None, compression='gzip',
                     usecols=[0, 3, 7], names=["Position", "Context", "Coverage"])

    # print a message indicating the file has been read
    print(f"File {sample_name} read")

    # filter for CpG context only
    df = df[df["Context"] == "CG"]

    # add column for the sample name
    df["Sample"] = sample_name

    # Append each sample's data to all_data, which is list
    subset_data.append(df)


# make sure len=number of samples read in - for retained samples
len(subset_data)
#165

# concatenate into giant dataframe (so that all_data is no longer an indexable list)
subset_data = pd.concat(subset_data)
#type(subset_data)

# when 'all_data' is a giant pandas DataFrame, index by column name
avg_coverage = subset_data["Coverage"].mean()

print(avg_coverage)
#20.846507652433736

# find SD
sd = subset_data.Coverage.std()
print(sd)
# 10.881626714867444

# find 95th percentile
high_cov = all_data.Coverage.quantile(q=0.95)
print(high_cov)
#51.0

# find 99th percentile
highest_cov = all_data.Coverage.quantile(q=0.99)
print(highest_cov)

#81.0

# Merge all data into a single DataFrame and change indexing of rows
coverage_data = pd.concat(all_data, ignore_index=True)

# Filter data to exclude high coverage cutoff
filtered_cov = all_data[all_data["Coverage"]<=88]

filtered_cov["Coverage"].mean()
#23.602394351882534

# Plot the first 5 samples
subset = all_data["Sample"].unique()[:5]

# Filter the DataFrame to include only the first 5 samples
subset_data = all_data[all_data["Sample"].isin(siub)]

# Plot using seaborn
plt.figure(figsize=(15, 8))
sns.violinplot(data=coverage_data[0'5, x="Sample", y="Coverage: nucleotides observed at CpGs", inner="quartile")
plt.xticks(rotation=90)
plt.xlabel("Sample")
plt.ylabel("Coverage at CG sites)"
#plt.ylim(0,150)
plt.title("Coverage Distribution Over Samples")
plt.tight_layout()

# Save the plot as a PNG file
# plt.savefig("coverage_violin_plot.png")
# plt.show()
