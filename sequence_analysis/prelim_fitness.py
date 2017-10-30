# imports
#############################################################
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
import plotly.plotly as py
import plotly.graph_objs as go

# load and manipulate dataframe
#############################################################

# load in dataframe of relative fitness values
d_rotation1 = pd.read_csv("jared_relative_fitness.csv")

# transpose - useful for some plotting applications
d_rotation2 = d_rotation1.transpose()

# turn -10 NaNs to actual NaNs for seaborn to mask in heatmap
# create seaborn mask
d_nan = d_rotation1.replace(-10, np.NaN)
mask1 = d_nan.isnull()

# just the c-terminus
d_nan_c_term = d_rotation1.drop(d_nan.columns[range(0,90)], axis=1).replace(-10, np.NaN)
mask2 = d_nan_c_term.isnull()

# collapsing dataframe
#############################################################

# find mean fitness for all values in the protein
mean_fit_all = []
for column in d_nan:
    mean_fit_all.append(d_nan[column].mean())

# maximum fitness for each residue
max_fit = []
for column in d_nan:
    max_fit.append(d_nan[column].max())

# min fitness for each residue
min_fit = []
for column in d_nan:
    min_fit.append(d_nan[column].min())

# plotting
#############################################################

# loads in wt sequence and changes to amino acid seq
wt = SeqIO.read("asyn_seq.txt", "fasta")
wt_aa = ((wt.seq.translate()))

# gets just the indices for charged amino acids
indices_neg = [i for i, x in enumerate(wt_aa) if x == "E" or x == "D"]
indices_pos = [i for i, x in enumerate(wt_aa) if x == "R" or x == "K"]

# subsets to just the fitness scores for aa's that are originally and saves as csv
pos_neg = d_rotation1[["E", "D", "R", "K"]]
pd.pos_neg.to_csv("pos_neg.csv")

# finds the neg residues and finds instances where neg residues are changed to positive
neg = d_rotation1.loc[indices_neg]
neg_fitness_to_pos = neg[["R", "K"]].replace(-10, np.NaN)

# finds positive residues and gets fitness of instances of positive changed to negative
pos = d_rotation1.loc[indices_pos]
pos_fitness_to_neg = pos[["E", "E"]].replace(-10, np.NaN)

# sends to plotly to make a heatmap of negative values changed to positive
data = [go.Heatmap( z=neg_fitness_to_pos.values.tolist(), colorscale='Viridis')]
py.iplot(data, filename='neg to pos')

# makes plotly heatmap of fitness pos to negative
data = [go.Heatmap(z=pos_fitness_to_neg.values.tolist(), colorscale='Viridis')]
py.iplot(data, filename='pos to neg')

# plots mean and max fitness for each residue
plt.plot(list(range(140)), max_fit)
plt.plot(list(range(140)), mean_fit_all)
plt.plot(list(range(140)), min_fit)
plt.show()

# dark style heatmap of all fitness values
sns.set_style("dark")
ax = sns.heatmap(d_rotation1, vmin=-2, vmax=2, mask=mask1, cmap=sns.diverging_palette(250, 15, s=75, l=40, center="dark", as_cmap=True), linewidths=.5)
ax = sns.heatmap(d_rotation1, vmin=-2, vmax=2, mask=mask1, cmap="hot", linewidths=.5)
plt.xlabel("residue number")
plt.ylabel("mutation")
plt.show()

