import pandas as pd
from os.path import join as opj
from fg_shared import _git, _fg_data

"""Prep the modules in a CSV from the supplementary info"""
project_folder = opj(_fg_data, 'SEATRAC', 'TB_hackday_2024')
liu_fn = opj(project_folder, 'seatrac-hackday-2024', 'liu_etal', 'Liu et al supplement', '1-s2.0-S266637912300215X-mmc2.xlsx')
module_fn = opj(project_folder, 'processed_data', 'liu_etal_modules.csv')
liu_df = pd.read_excel(liu_fn, skiprows=2).fillna('M8')

"""Merge with NHP data from Liu et al."""
mm2hum = pd.read_csv(opj(project_folder, 'processed_data', 'mmulatta_to_human.tsv'), sep='\t')
liu_df = pd.merge(liu_df, mm2hum[['Gene name', 'Human gene name']], how='left', on='Gene name')
liu_df.rename(columns={'Gene name':'gene', 'moduleID':'module', 'Gene stable ID':'geneid'})[['geneid', 'gene', 'module']].to_csv(module_fn, index=False)

"""Load required packages"""
import pandas as pd
import numpy as np
from os.path import join as opj
import sys
import seaborn as sns
from scipy import stats
import statsmodels.api as sm

"""Custom package used to define relative paths"""
from fg_shared import _git, _fg_data

project_folder = opj(_fg_data, 'SEATRAC', 'TB_hackday_2024')

"""Define file paths of the data we need to load"""
meta_fn =  opj(project_folder, 'processed_data', 'liu_etal_metadata.csv')
module_fn = opj(project_folder, 'processed_data', 'liu_etal_modules.csv')
cts_fn = opj(project_folder, 'processed_data', 'liu_etal_counts_long.csv')
immune_fn = opj(project_folder, 'processed_data', 'darrah_dose_immune_markers.csv')

"""Load the data and light processing"""
meta =  pd.read_csv(meta_fn)
modules = pd.read_csv(module_fn)
ncts = pd.read_csv(cts_fn)

imm = pd.read_csv(immune_fn)
imm = imm.rename(columns={'subjid':'animalid'})

"""Merge module names with the gene count data so that module scores can be computed for each module"""
ncts = pd.merge(ncts, modules, how='left', on='gene')

"""Module score is the average expression of all the genes in the module (using the log2-counts provided)"""
mod_scores = ncts.groupby(['sampleid', 'module'])['count'].agg('mean').reset_index()
mod_scores = pd.merge(mod_scores, meta, how='left', on='sampleid')


"""Quick plot of module scores by visit pre/post-BCG"""
sns.boxplot(data=mod_scores,
            x='module',
            y='count',
            hue='visit', 
            hue_order=['pre', 'd2', 'wk2', 'wk4', 'wk12'])

"""Estimate how M1 module scores at day 2 correlate with various immune responses measured later"""
mod = 'M1'
mod_visit = 'd2' #options ['pre', 'd2', 'wk2', 'wk12', 'wk4']
# imm_visit = 4 # wk4 but this is encoded in the key
key = 'Lg+ Marginal/CD4/IFNg Abs 4'

results = []
for key in imm.key.unique():
    """Join the module scores with the immune response dataset and look for correlations of M1 at d2 with a IV BCG response"""
    test_df = pd.merge(mod_scores.loc[(mod_scores['module'] == mod) & (mod_scores['visit'] == mod_visit)],
                      imm.loc[imm['key'] == key],
                      how='inner',
                      on='animalid')
    test_df = test_df.dropna(subset=['value', 'count'])
    rho, pvalue = stats.spearmanr(test_df['count'], test_df['value'])
    tmp = dict(module='M1',
               mod_visit=mod_visit,
               imm_key=key,
               n=test_df.shape[0],
               rho=rho,
               pvalue=pvalue)
    results.append(tmp)
results = pd.DataFrame(results).sort_values(by='pvalue')
results = results.assign(fdrq=sm.stats.multipletests(results['pvalue'].fillna(1), method='fdr_bh')[1])

results.head(20)

"""Example correlation scatter plot"""
mod = 'M1'
mod_visit = 'd2'
imm_key = 'PBMC/Mtb300 Marginal/CD4/CD154 12'


plot_df = pd.merge(mod_scores.loc[(mod_scores['module'] == mod) & (mod_scores['visit'] == mod_visit)],
                      imm.loc[imm['key'] == imm_key],
                      how='inner',
                      on='animalid')
plot_df = plot_df.dropna(subset=['value', 'count'])

sns.scatterplot(data=plot_df, x='count', y='value', hue='protect_outcome')
plt.xlabel(f'{mod} score at {mod_visit} visit')
plt.ylabel(imm_key)

