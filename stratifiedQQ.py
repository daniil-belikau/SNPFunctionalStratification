#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import assocplots.qqplot as qqplot
import numpy as np
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')

mpl.rcParams['figure.dpi']=100
mpl.rcParams['savefig.dpi']=100
mpl.rcParams['figure.figsize']=5.375, 5.375


# In[110]:


def plot_qq(data, title):
    """
    data: 2-tuple (list of p-values, list of labels)
    title: plot title
    
    return: void
    """
    
    qqplot.qqplot(data[0], 
           data[1], 
           color=['b','r', 'g'], 
           fill_dens=[0.2,0.2,0.2], 
           error_type='theoretical', 
           distribution='beta',
           title=title)

    mpl.pyplot.legend(loc=0)

    plt.savefig(title, dpi=300)
    plt.clf()


# In[33]:


def import_data(gwas, proms_path, ftypes, cells, states, chrom, bp, delim="\t", header=None):
    """
    gwas: 2-tuple (gwas path, important columns),
    proms_path: path to folder containing promoter data
    ftypes: list of data file types
    cells: list of cell lines
    states: list of possible epigenetic states
    delim: df sep
    header: header encoding in data files
    
    return: gwas dataframe, dictionary of promoter data
    """
    df_gwas = pd.read_csv(gwas[0], sep=delim)
    df_gwas = df_gwas[gwas[1]]
    
    proms = dict()
    for ftype in ftypes:
        proms[ftype] = dict()
        for cell in cells:
            proms[ftype][cell] = dict()
            for state in states:
                file_name = "_".join(["CapStarrseq", state, f"{cell}{'' if ftype == 'Proms' else '_Gwas_PromInfo'}.tsv"])
                file_path = os.path.join(proms_path, ftype, cell, state, file_name)
                df = pd.read_csv(file_path, sep=delim, header=header)
                df.sort_values(by=[chrom, bp], inplace=True)
                proms[ftype][cell][state] = df            
    print("Imported data")
    
    return df_gwas, proms


# In[4]:


def prepare_proms(proms, ftypes, cells, states):
    """
    proms: dict containing proms data
    ftypes: list of data file types
    cells: list of cell lines
    states: list of possible epigenetic states
    
    return: void
    """
    for ftype in ftypes:
        for cell in cells:
            for state in states:
                df = proms[ftype][cell][state]
                df[0] = df[0].apply(lambda s: s[3:])
                for col in range(3):
                    df[col] = df[col].apply(lambda s: int(s))
                if ftype == "vars":
                    df[3] = df[3].apply(lambda s: s[3:])
                    for col in range(3,6):
                        df[col] = df[col].apply(lambda s: int(s))
                    df[6] = df[6].apply(lambda s: s[:-2])
    print("Prepared proms")
    


# In[6]:


def add_to_gwas_dict(val, gwas_dict, chrom, bp):
    gwas_dict[val[chrom]].append(val[bp])
    

def create_gwas_dict(gwas, chrom, bp):
    """
    gwas: gwas dataframe
    chrom: name of chromosome column
    bp: name of base pair position column
    
    return: gwas dictionary
    """
    gwas.sort_values(by=[chrom, bp], inplace=True)
    gwas_dict = dict()
    for chromosome in gwas[chrom].unique():
        gwas_dict[chromosome] = list()
    gwas.apply(add_to_gwas_dict, args=(gwas_dict, chrom, bp), axis=1)
    print("Created GWAS dict")
    return gwas_dict
    


# In[34]:


def add_to_proms_dict(val, proms_dict, chrom, start, end):
    proms_dict[val[chrom]] += [val[start], val[end]]
    

def create_proms_dict(proms, chrom, start, end):
    """
    proms: promoter dataframe
    chrom: name of chromosome column
    start: name of start position column
    end: name of end position column
    
    return: proms dictionary
    """
    proms_dict = dict()
    for chromosome in proms[chrom].unique():
        proms_dict[chromosome] = list()
    proms.apply(add_to_proms_dict, args=(proms_dict, chrom, start, end), axis=1)
    print("Created prom dict")
    return proms_dict
    


# In[8]:


# def add_to_rsid_dict(val, rsid_dict, chrom, rsid):
#     rsid_dict[val[chrom]].append(val[rsid])
    

# def create_rsid_dict(proms, chrom, rsid):
#     """
#     proms: promoter dataframe
#     chrom: name of chromosome column
#     rsid: name of rsid column
    
#     return: rsid dictionary
#     """
#     proms.sort_values(by=[chrom, rsid], inplace=True)
#     rsid_dict = dict()
#     for chromosome in proms[chrom].unique():
#         rsid_dict[chromosome] = list()
#     proms.apply(add_to_rsid_dict, args=(rsid_dict, chrom, rsid), axis=1)
#     print("Created rsid dict")
#     return rsid_dict


# In[77]:


def get_prom_snps(gwas, gwas_dict, proms_dict, chrom_col, bp_col):
    """
    gwas: gwas dataframe
    gwas_dict: dictionary containing gwas data
    proms_dict: dictionary containing promoter data
    chrom_col: name of chromosome column
    bp_col: name of position column
    
    return: promoter snps dataframe
    """
    snps = list()
    for chrom in gwas[chrom_col].unique():
        if chrom in gwas_dict and chrom in proms_dict:
            print(chrom)
            for bp in gwas_dict[chrom]:
                groups = zip(*[iter(proms_dict[chrom])] * 2)
                if any(start <= bp <= end for start,end in groups):
                    snps.append((chrom, bp))  
                    
    snps_two = list()
    for snp in snps:
        snps_two.append(gwas[(gwas[chrom_col]==snp[0]) & (gwas[bp_col]==snp[1])])
    if len(snps_two) > 0:
        snps_two = pd.concat(snps_two)
    else:
        snps_two = [pd.DataFrame()]
    print("Collected prom snps based on bp")
    
    return snps_two


# In[10]:


def get_rsid_snps(gwas, prom_df, rsid_col_gwas, rsid_col_prom):
    """
    gwas: gwas dataframe
    prom_df: promoter dataframe 
    rsid_col_gwas: name of rsid column in gwas
    rsid_col_prom: name of rsid column in prom_df
    
    return: promoter snps dataframe
    """
    rsids = list(prom_df[rsid_col_prom].unique())
    snps = gwas[gwas[rsid_col_gwas].isin(rsids)]
    print("Collected prom snps based on rsid")
    
    return snps


# In[89]:


def distribute_snps(gwas, active, inactive, rsid_col):
    """
    gwas: gwas dataframe
    active: active epromoter dataframe 
    inactive: inactive epromoter dataframe
    rsid_col: name of rsid column
    
    return: non promoter snps, active snps, inactive snps dataframes
    """
    df_input = list()
    for temp in active:
        print(type(temp))
        if isinstance(temp, pd.core.frame.DataFrame):
            df_input.append(temp)
    active = pd.concat(df_input).drop_duplicates(rsid_col).reset_index(drop=True)
    df_input = list()
    for temp in inactive:
        print(type(temp))
        if isinstance(temp, pd.core.frame.DataFrame):
            df_input.append(temp)
    inactive = pd.concat(df_input).drop_duplicates(rsid_col).reset_index(drop=True)
    print("Combined all snp dfs")
    proms_all_snps = active.append(inactive)
    non_prom = gwas[~gwas[rsid_col].isin(proms_all_snps[rsid_col])]
    print("Separated non-prom snps")
    
    return non_prom, active, inactive


# In[87]:


def analyze_gwas(gwas_path, proms_path, ftypes, cells, states, chrom, bp, rsid, p, trait):
    gwas = (gwas_path, [chrom, bp, rsid, p])
    gwas, proms = import_data(gwas, proms_path, ftypes, cells, states, 0, 1)
    gwas_dict = create_gwas_dict(gwas, chrom, bp)
    prepare_proms(proms, ftypes, cells, states)
    for cell in cells:
        snps = dict()
        for state in states:
            snps[state] = list()
            for ftype in ftypes:
                if ftype == "Proms":
                    proms_dict = create_proms_dict(proms[ftype][cell][state], 0, 1, 2)
                    snps[state].append(get_prom_snps(gwas, gwas_dict, proms_dict, chrom, bp))
                else:
                    eproms_dict = create_proms_dict(proms[ftype][cell][state], 0, 1, 2)
                    gwas_proms_dict = create_proms_dict(proms[ftype][cell][state], 3, 4, 5)
                    snps[state].append(get_prom_snps(gwas, gwas_dict, eproms_dict, chrom, bp))
                    snps[state].append(get_prom_snps(gwas, gwas_dict, gwas_proms_dict, chrom, bp))
                    snps[state].append(get_rsid_snps(gwas, proms[ftype][cell][state], rsid, 6))
        non_prom, active, inactive = distribute_snps(gwas, snps["Active"], snps["Inactive"], rsid)
        data = ([non_prom[p], active[p], inactive[p]], ["Non-promoter", "Active", "Inactive"])
        plot_qq(data, f"{trait} {cell} Stratified QQ Plot")
    


# In[111]:


gwas = [("../data/gwas/Asthma2/Eur_clean_Asthma_data.txt", "chr", "position", "rsid", "European_ancestry_pval_rand", "Asthma European Ancestry"), 
        ("../data/gwas/Asthma2/Multi_clean_Asthma_data.txt", "chr", "position", "rsid", "Multiancestry_pval_rand", "Asthma Multiancestry"),
        ("../data/gwas/CD/clean_CD_data.txt", "CHR", "BP", "SNP", "P", "CD")]
proms_path = "../data/eprom_classification"
ftypes = ["Proms", "Vars"]
cells = ["HELA", "K562"]
states = ["Active", "Inactive"]

for g in gwas:
    analyze_gwas(g[0], proms_path, ftypes, cells, states, g[1], g[2], g[3], g[4], g[5])

