import pandas as pd
from os.path import join as opj
import numpy as np
import statsmodels.formula.api as smf
from scipy import stats
import sys

"""TODO:

 - Figure out what to put in the R script
 - Most promising is probably to find gran-associated CD4 T cell genes
 - Then test which are associated with gran CFU
 - Validation can be with Bromley pseudo bulk?
 - Hack could be to find gene sets instead to boost signal
 """

_git = r'C:\Andrew\gitrepo'
data_folder = opj(_git, 'seatrac-hackday-2024', 'foreman_etal')

tpm = pd.read_csv(opj(data_folder, 'ge_tpm.tsv'), sep='\t')
edds = pd.read_csv(opj(data_folder, 'edds.tsv'), sep='\t')

keep_cols = ['ColumnName', 'biosample', 'biosample_name', 'subject', 'compartment', 'sort', 'gran',
               'Granuloma__', 'CD4_T_cells', 'CD8_T_cells', 'CFU',
               'CD11b_Mtb_', 'CD11b_Mtb_0', 'condition']

edds = edds[keep_cols].rename(columns={'Granuloma__':'Granuloma_num',
                                        'ColumnName':'sampleid'})
edds = edds.assign(logCFU=np.log10(edds['CFU'] + 1))

ygenes = ['DDX3Y', 'UTY', 'ZFY']
ytpm = tpm.set_index('gene_id').loc[ygenes, edds['sampleid']]
sex = np.log10(ytpm.sum(axis=0) + 1)
sex.name = 'sex'
sex = sex.map(lambda v: 'M' if v < 0.5 else 'F')
edds = edds.set_index('sampleid').join(sex).reset_index()

edds.to_csv(opj(data_folder, 'foreman_etal_meta.csv'), index=False)

tpm = tpm.loc[:, ['gene_id'] + list(edds['sampleid'])]
tpm.loc[:, edds['sampleid']] = np.log2(tpm.loc[:, edds['sampleid']] + 1)
tpm.to_csv(opj(data_folder, 'foreman_etal_counts.csv'), index=False)

"""Look for genes that are associated with granulomas
SIgnificant finding that there are granuloma associated CD4 T cell genes"""
meta = edds.loc[edds['condition'].isin(['CD4_gran', 'CD4_PBMC'])].set_index('sampleid')
ltpm = tpm.set_index('gene_id')
ltpm = ltpm.loc[(ltpm > 1).mean(axis=1) > 0.6]
genes = ltpm.index
mdf = meta.join(ltpm.T, how='left')

ind = mdf['condition'] == 'CD4_gran'

res = []
for g in genes:
    gran_expr = mdf.loc[ind, g]
    pbmc_expr = mdf.loc[~ind, g]
    stat, pvalue = stats.mannwhitneyu(gran_expr, pbmc_expr)
    res.append({'gene':g,
        'stat':stat,
        'pvalue':pvalue,
        'assoc':'GRAN' if np.mean(gran_expr) > np.mean(pbmc_expr) else 'PBMC'})
res = pd.DataFrame(res).sort_values(by='pvalue')
res = res.assign(FDRq=sm.stats.multipletests(pvals=res['pvalue'], method='fdr_bh')[1])
print(res.head(20))

"""Test out correlation of genes with granuloma CFU
Maybe the exercise would be to redo this analysis using modeling instead of rank-corr?"""
cd4_genes = ['KLRB1', 'CD40LG', 'S100A11', 'S100A4', 'IL26', 'BATF']
cd8_genes = ['APOBEC3G', 'IFNG', 'TNF', 'CCL1','CCL20']
sig_genes = res.loc[res['FDRq'] < 0.1, 'gene'].tolist()

meta = edds.loc[edds['condition'].isin(['CD4_gran'])].set_index('sampleid')
ltpm = tpm.set_index('gene_id')
mdf = meta.join(ltpm.loc[sig_genes].T, how='left')

cfu_res = []
for gene in sig_genes:
    """OLS model"""
    model = smf.ols(formula=f'logCFU ~ {gene}',
                        data=mdf)

    summ = model.fit()
    lb, ub = summ.conf_int().loc[gene]
    rho, corr_pvalue = stats.spearmanr(mdf['logCFU'], mdf[gene])
    cfu_res.append({'gene':gene,
                    'est':summ.params[gene],
                    'est_lb':lb,
                    'est_ub':ub,
                    'pvalue':summ.pvalues[gene],
                    'corr_pvalue':corr_pvalue,
                    'rho':rho})

cfu_res = pd.DataFrame(cfu_res).sort_values(by='pvalue')
cfu_res = cfu_res.assign(FDRq=sm.stats.multipletests(pvals=cfu_res['pvalue'], method='fdr_bh')[1])
print(cfu_res.head(20))


"""Mixed-effects model (not much variance modeled within animal since grans are independent)"""
model = smf.mixedlm(formula=f'logCFU ~ {gene} + condition',
                    data=mdf,
                    groups='subject',
                    re_formula='1')

res = model.fit()
print(res.summary())

sns.jointplot(data=mdf, x='logCFU', y=gene)
print(stats.spearmanr(mdf['logCFU'], mdf[gene]))