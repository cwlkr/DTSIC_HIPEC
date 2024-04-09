import os
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import patsy

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import logrank_test, pairwise_logrank_test, multivariate_logrank_test
import seaborn as sns

from scipy.stats import zscore, pearsonr
from sklearn.model_selection import LeaveOneOut, StratifiedKFold, KFold
from scipy.spatial.distance import cosine

survevent = 'status'
survtime = 'time'
data_org = pd.read_csv('dwls_results.csv')  # should include pfs time and status

## immune only, original cluster names from deconv
features_selected = [
'Immune_cells_1_makrophages',
 'Immune_cells_2_cd8t_cell',
 'Immune_cells_3_b_cells',
 'Immune_cells_4_makrophages',
 'Immune_cells_5_b_cells']

data_org = data_org[data_org['Moment'] == 'resection']  # 149
data_org = np.where(data_org['treatment'] == 'CRS', 0, 1) # binarise for easier CPH calculations; is equiavalent to categorical case

# mean (the two) patients with multiple resection samples
agg_f = {**{i:'mean' for i in features_selected},**{i: lambda x: x.iloc[0] for i in data_org.columns if i not in features_selected}}
data_org = data_org.groupby('PtNr').agg(agg_f).reset_index(drop=True)
assert not np.any(data_org['PtNr'].duplicated()), 'duplicated patients'

formula = f"0+ status + time + ({' + '.join([f for f in features_selected])})*treatment"

# data_org[features_selected] = data_org[features_selected].apply(lambda x:zscore(x)) # optinal zscore; CPH is scaling invariant
outdir = 'figure_output/'

data_ = patsy.dmatrix(formula, data_org, return_type="dataframe")
data_ = data_.reset_index(drop=True)

features_selected = data_.columns[np.where((data_.columns != 'time') & (data_.columns != 'status'))].to_list()

idx_t = np.asarray(features_selected)[np.argwhere(np.array(features_selected) == 'treatment')].squeeze().tolist()
idx_vars = np.asarray(features_selected)[np.argwhere(np.array([False  if ':treatment' in i or i == 'treatment' else True for i in features_selected]))].squeeze().tolist()
idx_int = np.asarray(features_selected)[np.argwhere(np.array([True  if ':treatment' in i else False for i in features_selected]))].squeeze().tolist()


# n_splits = 25
# cv =KFold(n_splits=n_splits, shuffle=True).split(data_)
model_params = []
dts_scores = np.zeros((len(data_)))

cv = LeaveOneOut().split(data_)
for i, (train_idx, test_idx) in enumerate(cv):
    X_train, X_test = data_.loc[train_idx][features_selected + ['time', 'status']], data_.loc[test_idx][features_selected + ['time', 'status']]
    assert len(X_train.merge(X_test)) == 0  # check that there is no overlap!
    cph = CoxPHFitter()
    cph.fit(X_train, duration_col=survtime, event_col=survevent)
    coef_t, ct_HR, coef_int = cph.params_[idx_t], X_test[idx_vars], cph.params_[idx_int]
    dts_scores[test_idx] = (-( coef_t + (np.asarray(ct_HR) @ np.asarray(coef_int)[None,:].T))).squeeze()
    model_params.append(np.asarray(cph.params_))


coefs = pd.DataFrame(np.array(model_params), columns=cph.params_.keys())

mapping = {'Immune_cells_1_makrophages:treatment': 'Macrophages SPP1+:T',
 'Immune_cells_4_makrophages:treatment': 'Macrophages:T',
 'Immune_cells_3_b_cells:treatment': 'B-cells SSR4+:T',
 'Immune_cells_2_cd8t_cell:treatment': 'T-cells:T',
 'Immune_cells_5_b_cells:treatment': 'B-cells:T'}
# print(coefs.mean().sort_values())
mean_c = coefs.mean().sort_values()
coef_means = coefs.reindex(coefs.mean().sort_values().index, axis=1)
coefs = coef_means.melt()
coef_p = coefs.loc[coefs.covariate.str.contains(':'),].replace(mapping)

font = {'size': 17}
matplotlib.rcdefaults()
matplotlib.rc('font', **font)
fig, ax = plt.subplots(figsize=(9,5))
c_d = coef_p.rename(columns={'covariate': 'Cell type', 'value':'log(HR)'})
c_d['Cell type'] = c_d['Cell type'].astype("string")
ax = sns.barplot(c_d, y= 'Cell type', x = 'log(HR)', color = 'grey', capsize=0.1)
ax.set_ylabel('')
# ax.set_title('Interaction coefficients from leave-one-out validation')
data_.apply('max', axis=0)
ax.axvline(0, color=".3", dashes=(2, 2))
plt.tight_layout()
# fig.savefig(f'{outdir}/coef_plot.png', dpi=400)
plt.show()
matplotlib.rcdefaults()

def plot_kaplan_c(f, surv_frame, strat, t=None, ll = ['CRS', 'HIPEC']):
    surv_frame = surv_frame[[f, survtime, survevent,strat]].dropna()
    fig, ax = plt.subplots(1, figsize=(10,10))
    ms = []
    for i in np.unique(surv_frame[strat]):
        surv_frame_s = surv_frame.loc[surv_frame[strat] == i,  ['time', 'status', f, strat]]
        if t is None:  # t per group!!!!
            t_ = np.median(surv_frame_s[f])
        else:
            t_ = t
        fiture =  (surv_frame_s[f] >= t_)
        T,E = surv_frame_s[[survtime]]/365, surv_frame_s[[survevent]]
        kmf_high = KaplanMeierFitter()
        ax = kmf_high.fit(T[fiture], event_observed=E[fiture], label=f'{ll[int(i)]}_{f}_high').plot_survival_function(ax=ax, show_censors=True)
        ms.append(kmf_high)
        kmf_low = KaplanMeierFitter()
        ax = kmf_low.fit(T[~fiture], event_observed=E[~fiture], label=f'{ll[int(i)]}_{f}_low').plot_survival_function(ax=ax, show_censors=True)
        ms.append(kmf_low)
        print(logrank_test(T[fiture],T[~fiture], E[fiture],E[~fiture]), f':treat {ll[int(i)]}')
    add_at_risk_counts(*ms, ax=ax)
    ax.legend(loc='upper center', ncols=4)
    ax.set_xlabel('Time (Years)')
    ax.set_ylabel('Survival probability')
    # plt.title(title)
    plt.tight_layout()
    plt.show()
    return fig

strat = 'treatment'
data_['DTS'] = dts_scores

plot_kaplan_c('DTS', data_, strat='treatment')

stats_frame = data_.copy()
stats_frame['DTS'] = zscore(stats_frame['DTS'])
stats_frame_ = patsy.dmatrix(f"0+ status + time + DTS * treatment", stats_frame, return_type="dataframe")
stats_frame_ = stats_frame_.reset_index(drop=True)

cph_r = CoxPHFitter()
cph_r.fit(stats_frame_[['DTS', 'treatment', 'DTS:treatment', survevent, survtime]], duration_col=survtime, event_col=survevent)
cph_r.print_summary()### feature selection done on all cases using cross validation.
cph_r.plot(hazard_ratios=True)
plt.show()
print(cph_r.log_likelihood_ratio_test())



a, b = np.median(stats_frame_[stats_frame_['treatment'] == 0]['DTS']),  np.median(stats_frame_[stats_frame_['treatment'] == 1]['DTS'])
stats_frame_['groups'] = np.where((stats_frame_['treatment'] == 1) & (stats_frame_['DTS'] > b) , 'HIPEC_DTS_high',
                         np.where((stats_frame_['treatment'] == 1) & (stats_frame_['DTS'] <= b) , 'HIPEC_DTS_low',
                         np.where((stats_frame_['treatment'] == 0) & (stats_frame_['DTS'] < a) , 'CRS_DTS_low', 'CRS_DTS_high')))

font = {'size': 15}
matplotlib.rcdefaults()
matplotlib.rc('font', **font)

def plot_kaplan_cat(f, surv_frame, survevent='status', survtime='time'):
    feature = f
    surv_frame = surv_frame[[f, survtime, survevent]].dropna()

    fig, ax = plt.subplots(1, figsize=(10,10))
    cats = np.unique(surv_frame[f])

    l = []
    for i in np.roll(cats, 1):
        fiture =  surv_frame[f] == i
        T = surv_frame[[survtime]]/365
        E = surv_frame[[survevent]]
        kmf_high = KaplanMeierFitter()
        ax = kmf_high.fit(T[fiture], event_observed=E[fiture], label=i.split('_')[0]).plot_survival_function(ax=ax, show_censors=True)
        l.append(kmf_high)

    add_at_risk_counts(*l, ax=ax)
    handles, labels = ax.get_legend_handles_labels()
    order = [0,1]
    ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc='upper center', ncols=2)
    ax.set_xlabel('Time (Years)')
    ax.set_ylabel('Progression-free survival (probability)')
    plt.tight_layout()
    plt.show()
    result = multivariate_logrank_test(surv_frame[survtime], surv_frame[f], surv_frame[survevent])
    result.print_summary()
    return fig

stats_frame_['groups'] = stats_frame_['groups'].replace({'HIPEC_DTS_high':'CRS-HIPEC_DTS_high',
                                                'HIPEC_DTS_low':'CRS-HIPEC_DTS_low'})

f1 = plot_kaplan_cat('groups', stats_frame_.loc[stats_frame_['groups'].str.contains('igh')])
f2 = plot_kaplan_cat('groups', stats_frame_.loc[stats_frame_['groups'].str.contains('low')])

# plot_kaplan_c('DTS', stats_frame_.loc[stats_frame_['treatment']==1,:], strat, t=None)
# plot_kaplan_c('DTS', stats_frame_.loc[stats_frame_['treatment']==0,:], strat, t=None)
f1.savefig(os.path.join(outdir, f'KM_DTS_IC_HIGH.png'), dpi=400)
f2.savefig(os.path.join(outdir, f'KM_DTS_IC_LOW.png'), dpi=400)

print(pairwise_logrank_test(stats_frame_[survtime], stats_frame_['groups'], stats_frame_[survevent]))

data_['DTS_m'] = np.where(data_['DTS'] < np.median(data_['DTS']), 'LOW' , 'HIGH' )
data_ = data_.reset_index()

plot_df = pd.melt(data_, value_vars= [f for f in features_selected if 'tr' not in f], id_vars=['index', 'DTS_m', 'treatment'], var_name = 'cell_type', value_name='abundance')

mapping = {'Immune_cells_1_makrophages': 'Macrophages SPP1+', 'Immune_cells_4_makrophages':'Macrophages',
            'Immune_cells_3_b_cells':'B-cells SSR4+', 'Immune_cells_2_cd8t_cell': 'T-cells', 'Immune_cells_5_b_cells':'B-cells'}

feature_ordered = [i.split(':')[0] for i in mean_c.keys() if ':' in i]
plot_df = pd.melt(data_[feature_ordered + ['treatment', 'DTS_m']], value_vars= feature_ordered, id_vars=[ 'DTS_m', 'treatment'], var_name = 'cell_type', value_name='abundance')

data_.columns
for k,v in mapping.items():
    data_[v] = data_[k]

plot_df.cell_type = plot_df.cell_type.apply(lambda x: mapping[x])


font = {'size': 20}
matplotlib.rcdefaults()
matplotlib.rc('font', **font)
fig, ax = plt.subplots(1, figsize=(10,10))
g = sns.boxplot(data = plot_df, x = 'cell_type', y='abundance', hue='DTS_m' , palette=[sns.color_palette("tab10")[3], sns.color_palette("tab10")[0]], ax=ax)
plt.xticks(rotation=45)
plt.title('')
plt.xlabel(' ')
plt.ylabel('Cell type abundance')
legend = ax.legend()
legend.texts[0].set_text(r"DTS$\mathregular{^{IC}}$-LOW")
legend.texts[1].set_text(r"DTS$\mathregular{^{IC}}$-HIGH")
plt.tight_layout()
fig.savefig(f'{outdir}/abundance_box_ordered_ic.png', dpi=400)
plt.show()


fig, ax = plt.subplots(1, figsize=(10,10))
g = sns.violinplot(data = plot_df, x = 'cell_type', y='abundance', hue='DTS_m' , palette=[sns.color_palette("tab10")[3], sns.color_palette("tab10")[0]], ax=ax)
plt.xticks(rotation=45)
plt.title('')
plt.xlabel(' ')
plt.ylabel('Cell type abundance')
legend = ax.legend()
legend.texts[0].set_text("DTS-LOW")
legend.texts[1].set_text("DTS-HIGH")
plt.tight_layout()
fig.savefig(f'{outdir}/abundance_violin_ic.png', dpi=400)

data_['DTS'] = zscore(data_['DTS'])
data_['time'] = data_['time']/365

# plot interaction
fig, ax = plt.subplots(figsize=(10,10))
sns.regplot(x="time", y="DTS", data=data_.loc[data_['status']==1].loc[data_['treatment']==1], ci=0, ax=ax, label="HIPEC", truncate=False, order=1)
sns.regplot(x="time", y="DTS", data=data_.loc[data_['status']==1].loc[data_['treatment']==0], ci=0, ax=ax, label="CRS", truncate=False, order=1)
ax.legend()
ax.set_xlabel('Time (Years)')
ax.set_ylim(-3,3)
plt.show()

# fig.savefig(f'{outdir}/scatterplot_time_DTS.png', dpi=400)

# other clean clusters are in clean_clusters.py file
meta_data_deseq = data_org.copy()
meta_data_deseq['DTS'] = dts_scores
meta_data_deseq['DTS_cat'] = np.where(meta_data_deseq['DTS'] < np.quantile(meta_data_deseq['DTS'], 0.25), 'LOW' , np.where(meta_data_deseq['DTS'] > np.quantile(meta_data_deseq['DTS'], 0.75), 'HIGH', 'MID'))
meta_data_deseq['DTS_bin'] = np.where(meta_data_deseq['DTS'] < np.quantile(meta_data_deseq['DTS'], 0.5),'LOW','HIGH')
# save dts with categories
# meta_data_deseq[['V1', 'DTS', 'DTS_cat', 'DTS_bin','treatment', 'status', 'time']].to_csv('/data/nki_ovarian_scs/deseq_dts_metadata.csv')
