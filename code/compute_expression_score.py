#!/usr/bin/python
import os
import sys
import logomaker
import pandas as pd
import numpy as np

def extract_Ab_info(filename):
  datasheet = pd.read_excel(filename, sheet_name='Sheet1')
  ab_dict = {}
  for index, info_dict in datasheet.iterrows():
    for i in info_dict.keys():
      info_dict[i] = str(info_dict[i])
    ab_dict[index]  = info_dict
  return ab_dict

def count_to_freq(df, colname):
    df[colname+'_freq'] = (df[colname]+1)/(df[colname].sum()+len(df))
    return (df)

def exp_score_calculate(df, rep, freq_cutoff):
    print (rep)
    df['Exp_weight_'+rep] = np.log10(df['PE_bin2_Expression_'+rep+'_freq']/df['PE_bin0_Expression_'+rep+'_freq'])
    df_high_freq = df[df['avg_input_freq'] >= freq_cutoff]
    print(len(df_high_freq))
    w_summary = df_high_freq.groupby('class')['Exp_weight_'+rep].mean()
    w_summary = w_summary.reset_index()
    w_silent   = (float(w_summary.loc[w_summary['class']=='WT']['Exp_weight_'+rep]))
    w_nonsense = (float(w_summary.loc[w_summary['class']=='nonsense']['Exp_weight_'+rep]))
    df['Exp_score_'+rep] = (df['Exp_weight_'+rep]-w_nonsense)/(w_silent-w_nonsense)
    print ('w_nonsense', w_nonsense)
    print ('w_WT', w_silent)
    return (df)

def classify_ab(ID, ab_name_flu, ab_geneID_flu, ab_geneID_other):
  if ID[-2::] == '_1':
    ID = ID[0:-2]
  if 'CR9114_WT' in ID:
    return 'WT'
  elif 'CR9114' in ID and ID[-1] == '_':
    return 'nonsense'
  elif 'CR9114' in ID:
    return 'CR9114 mutant'
  elif ID in ab_name_flu:
    return 'influenza'
  elif ID in ab_geneID_flu:
    return 'influenza'
  elif ID in ab_geneID_other:
    return 'others'
  else:
    print ('unclassified antibody: %s' % ID)
    return 'unclassified'

def count_to_express_score(outfile, file_count, ab_name_flu, ab_geneID_flu, ab_geneID_other):
  df = pd.read_csv(file_count)
  df = df.drop('Unnamed: 0', axis=1)
  df = df.fillna(0)
  df = df.loc[df['concentration']=='A']
  df['class'] = df['ID'].apply(classify_ab, args=(ab_name_flu, ab_geneID_flu, ab_geneID_other))
  df = df.loc[df['class']!='unclassified']
  df = count_to_freq(df, 'Input1_Expression_1')
  df = count_to_freq(df, 'Input2_Expression_2')
  df = count_to_freq(df, 'PE_bin0_Expression_1')
  df = count_to_freq(df, 'PE_bin1_Expression_1')
  df = count_to_freq(df, 'PE_bin2_Expression_1')
  df = count_to_freq(df, 'PE_bin0_Expression_2')
  df = count_to_freq(df, 'PE_bin1_Expression_2')
  df = count_to_freq(df, 'PE_bin2_Expression_2')
  df['avg_input_freq'] = (df['Input1_Expression_1_freq'] + df['Input2_Expression_2_freq']) /2
  df = exp_score_calculate(df, '1', 0)
  df = exp_score_calculate(df, '2', 0)
  df['Exp_score'] = (df['Exp_score_1'] + df['Exp_score_2'])/2
  print ('writing: %s' % outfile)
  df.to_csv(outfile, sep=",", index=False)

def main():
  outfile = 'result/CDRH3_express_table_summary.csv'
  file_count = 'result/CDRH3_count_table.csv'
  file_flu169 = 'doc/Flu_1-69.xlsx'
  file_all169 = 'doc/GenBank_1-69.xlsx'
  ab_dict_flu = extract_Ab_info(file_flu169)
  ab_dict_all = extract_Ab_info(file_all169)
  ab_name_flu = set([ab_dict_flu[i]['Name'] for i in ab_dict_flu])
  ab_geneID_flu = set([ab_dict_flu[i]['VH Genbank ID'] for i in ab_dict_flu])
  ab_geneID_all = set([ab_dict_all[i]['Genbank ID'] for i in ab_dict_all])
  ab_geneID_other = ab_geneID_all-ab_geneID_flu
  count_to_express_score(outfile, file_count, ab_name_flu, ab_geneID_flu, ab_geneID_other)

if __name__ == "__main__":
  main()
