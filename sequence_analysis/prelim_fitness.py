
import pickle
import pandas
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO


# load in dataframe
d_rotation1 = pickle.load(open("buds_fit_df1.pkl", "rb"))

# transpose - maybe useful?
d_rotation2 = d_rotation1.transpose()

# turn -10 NaNs to actual NaNs for seaborn to mask in heatmap
# create seaborn mask
d_nan = d_rotation1.replace(-10, np.NaN)
mask1 = d_nan.isnull()

# just the c-terminus
d_nan_c_term = d_rotation1.drop(d_nan.columns[range(0,90)], axis=1).replace(-10, np.NaN)
mask2 = d_nan_c_term.isnull()
#
res_of_interest = [120, 121, 122, 127, 128, 129, 134]

mean_fit_interesting = []
for res in res_of_interest:
    mean_fit_interesting.append(d_nan[res].mean())

mean_fit_all = []
for column in d_nan:
    mean_fit_all.append(d_nan[column].mean())

max_fit = []
for column in d_nan:
    max_fit.append(d_nan[column].max())

min_fit = []
for column in d_nan:
    min_fit.append(d_nan[column].min())

wt = SeqIO.read("asyn_seq.txt", "fasta")
wt_aa = list((wt.seq.translate()))

# print(wt_aa[90:])

indices_neg = [i for i, x in enumerate(wt_aa) if x == "E" or x == "D"]
indices_pos = [i for i, x in enumerate(wt_aa) if x == "R" or x == "K"]
print(indices_neg)



# print(np.argwhere(np.isnan(min_fit)))
# print(len(wt_list))
#
# plt.plot(res_of_interest, mean_fit_interesting, 'o')
# plt.plot(list(range(140)), max_fit)
# plt.plot(list(range(140)), mean_fit_all)
# plt.plot(list(range(140)), min_fit)
# plt.show()


sns.heatmap(d_rotation1, vmin=-8, vmax=2, mask=mask1, cmap="coolwarm")
plt.show()

