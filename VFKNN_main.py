import pandas as pd
import os, re, pickle, copy
import numpy as np
from sklearn import preprocessing
import warnings
warnings.filterwarnings('ignore')
# %%
filepath = 'E:/bioinformatics_linux/globalsoil/temp/readsvftpm/'
kraphlanpath = 'E:/bioinformatics_linux/globalsoil/temp/kraphlan/'
modelpath='E:/bioinformatics_linux/db/pathogen/NPDPsoilpathogendataset/'
out = 'E:/bioinformatics_linux/globalsoil/temp/knnpred/onMetaphlan/'
# %%
with open(modelpath + "knnmodel_trained.pickle", "rb") as file:
    clf, _, _ = pickle.load(file)
with open(out + "colnames.pickle", "rb") as file:
    colnames = pickle.load(file)
#
columns = ['qseqid', 'VFGID', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'qlen', 'sstart', 'send',
           'slen', 'evalue', 'bitscore', 'qcov', 'scov', 'nreads_within_s', 'averl_within_s', 'rpk_within_s',
           'tpm_within_s', 'taxid', 'nReads', 'n16S', 'nCell', 'VFID', 'VFCID.x', 'VFCID.y', 'proteinid', 'VF_Name',
           'VFcategory', 'Function', 'Mechanism', 'Reference', 'VF_FullName', 'Bacteria', 'Characteristics',
           'Structure', 'sseqlength']

def pred_fuc(df, sample):
    grouped = df.groupby(by=['taxid', 'VFcategory', 'VF_Name', 'VFGID'])['tpm_within_s'].mean()
    df1 = pd.DataFrame(grouped).rename(columns={'0': 'tpm_within_s'}).reset_index()
    grouped = df1.groupby(by=['taxid', 'VFcategory', 'VF_Name'])['tpm_within_s'].sum()
    df1 = pd.DataFrame(grouped).rename(columns={'0': 'tpm_within_s'}).reset_index()
    r1 = copy.deepcopy(df1)
    GCFdata = r1.groupby(by=['taxid'])['tpm_within_s'].sum().reset_index().rename(columns={'tpm_within_s': 'all'})
    ch1 = GCFdata._append(pd.DataFrame(columns=r1["VFcategory"].unique()))._append(pd.DataFrame(columns=r1["VF_Name"].unique()))
    GCFcat = r1.groupby(by=['taxid', 'VFcategory'])['tpm_within_s'].sum().reset_index().rename(
        columns={'tpm_within_s': 'vfcattpm'})
    for i in range(0, len(GCFcat)):
        ch1.loc[ch1[ch1.taxid == GCFcat.loc[i, 'taxid']].index.tolist()[0], GCFcat.loc[i, 'VFcategory']] = GCFcat.loc[i, 'vfcattpm']
    for i in range(0, len(r1)):
        ch1.loc[ch1[ch1.taxid == r1.loc[i, 'taxid']].index.tolist()[0], r1.loc[i, 'VF_Name']] = r1.loc[i, 'tpm_within_s']
    ch2 = ch1.fillna(0).sort_values(['taxid'], ascending=[True])
    ch_name = pd.DataFrame(columns=colnames.to_list())
    common_columns = colnames.intersection(ch2.columns)
    for i in ch2.columns.to_list():
        if i not in common_columns:
            print(i)
    for i in colnames.to_list():
        if i in ch2.columns.to_list():
            ch_name[i] = ch2[i]
        else:
            ch_name[i] = 0
    ch2 = copy.deepcopy(ch_name)

    # prediction
    X = np.array(ch2.iloc[:, 1:])
    X = preprocessing.scale(X)
    y_pred = clf.predict(X)
    y_score = clf.predict_proba(X)

    result = pd.DataFrame(np.hstack((np.array(ch2.iloc[:, 0:1]), y_pred.reshape(-1, 1), y_score[:, 1].reshape(-1, 1)))).rename(
        {0: 'taxid', 1: 'type', 2: 'possibility'}, axis=1)
    result['taxid'] = result['taxid'].astype('int')
    result['type'] = result['type'].astype('int').astype('str')

    result.to_csv(out + sample + '_pre.tsv', sep='\t', header=True, index=False)
    ch2.to_csv(out + sample + '_pre_ch2.tsv', sep='\t', header=True, index=False)
    return result, ch2

def pre_tax_relab_df_fuc(filepath, f, sample, kraphlanpath, out, ab=''):

    fp = os.path.join(filepath, f)
    data = pd.read_table(fp, header=0, names=columns,usecols=['taxid', 'VFcategory', 'VF_Name', 'VFGID', 'tpm_within_s', 'nCell'])
    #
    data2 = pd.read_table('E:/bioinformatics_linux/globalsoil/temp/kraphlan/'+sample+'_tavg_g_profiled_metagenome.txt', header=4, names=['clade_name','NCBI_tax_id','relative_abundance','additional_species'],usecols=['clade_name', 'NCBI_tax_id'])
    data2=data2[data2['clade_name'].str.contains('\|t__')]
    data2[['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species','strain']]=data2['NCBI_tax_id'].str.split('|', expand=True)
    data=data[data['taxid'].astype(str).isin(data2['species'].tolist())]
    data = data[data['nCell'] != 0]
    if data['taxid'].nunique() < 2:
        return
    data['tpm_within_s'] = data['tpm_within_s'] / data['nCell']

    # predict
    if sample + '_pre.tsv' in os.listdir(out) and sample + '_pre_ch2.tsv' in os.listdir(out):
        predf = pd.read_table(out + sample + '_pre.tsv', header=0, dtype={'type': str})
        # ch2 = pd.read_table(out + sample + '_pre_ch2.tsv', header=0)
    else:
        predf, ch2 = pred_fuc(data, sample)

    # taxon
    if sample + '_pre_tax_relab' + ab + '.tsv' in os.listdir(out):
        tax_relab = pd.read_table(out + sample + '_pre_tax_relab' + ab + '.tsv', header=0)
    tax_relab = pd.read_table(kraphlanpath + sample + '_tax_relab' + ab + '.tsv', header=0)

    # merge
    pre_tax_df = pd.merge(tax_relab, predf, left_on='taxid', right_on='taxid', how='left')
    pre_tax_df['type'] = pre_tax_df['type'].fillna('0')
    pre_tax_df['sample'] = sample

    pre_tax_df.to_csv(out + sample + '_pre_tax_relab' + ab + '.tsv', sep='\t', header=True,
                      index=False)
    return

# %%

'''vfknn（tax,relab,prediction）'''
count = 0
for f in os.listdir(filepath):
    if f.endswith('_reads_vf_TPM.tsv'):
        sample = re.sub(r'_.+.tsv', '', str(f))
        count += 1
        print(sample, count)
        pre_tax_relab_df_fuc(filepath, f, sample, kraphlanpath, out, ab='_ab01')

