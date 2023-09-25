#!/usr/bin/python
import os
import sys
import glob
import logomaker
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def extract_Ab_info(filename):
  datasheet = pd.read_excel(filename, sheet_name='Sheet1')
  ab_dict = {}
  for index, info_dict in datasheet.iterrows():
    for i in info_dict.keys():
      info_dict[i] = str(info_dict[i])
    ab_dict[index]  = info_dict
  return ab_dict

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
    
def process_KD_file(file_KD, ab_name_flu, ab_geneID_flu, ab_geneID_other, epi_dict, seq_dict, outfile):
  df = pd.read_csv(file_KD)
  df['class'] = df['ID'].apply(classify_ab, args=(ab_name_flu, ab_geneID_flu, ab_geneID_other))
  df['epi'] = df['ID'].map(epi_dict)
  df['seq'] = df['ID'].map(seq_dict)
  df['CDRH3_len'] = df['seq'].str.len()
  df['resi95'] = df['seq'].str.slice(2,3)
  df['resi97'] = df['seq'].str.slice(4,5)
  df['resi98'] = df['seq'].str.slice(5,6)
  df.to_csv(outfile, index=False, na_rep='nan')
  print ("written: %s" % outfile)
  return df

def parse_epi_file(infile):
  epi_dict = {}
  infile = open(infile,'r')
  for line in infile.readlines():
    if 'ID' in line and 'epi' in line: continue
    else:
      line = line.rstrip().rsplit("\t")
      ID   = line[0]
      epi  = line[1]
      if epi == 'HA:Stem': epi  = 'HA stem'
      elif epi == 'HA:Unk': epi = 'HA unknown'
      else: 'Influenza (other)'
      epi_dict[ID] = epi
  return epi_dict

def seq_to_dict(file_seq):
  infile = open(file_seq, 'r')
  seq_dict = {}
  for line in infile.readlines():
    if 'Genbank ID' in line: continue
    else:
      line = line.rstrip().rsplit("\t")
      ID   = line[0]
      try:
        seq  = line[1]
      except:
        continue
      seq_dict[ID] = seq
      seq_dict[ID+'_1'] = seq
  return (seq_dict)

def make_sequence_logo(sequence_list, CDRH3_len, figname):
  sequence_list = [seq for seq in sequence_list if len(seq)==CDRH3_len]
  logo_width = (CDRH3_len-0)*0.6
  logo_width = 11 if logo_width > 11 else logo_width
  fig, ax = plt.subplots(1,1,figsize=[logo_width,2])
  seqlogo_matrix = logomaker.alignment_to_matrix(sequence_list)
  seqlogo = logomaker.Logo(seqlogo_matrix, font_name="Arial", color_scheme="weblogo_protein", width=1, ax=ax)
  seqlogo.style_spines(visible=False)
  seqlogo.ax.set_xticks([])
  seqlogo.ax.set_yticks([])
  seqlogo.ax.spines['bottom'].set_visible(False)
  seqlogo.fig.tight_layout()
  plt.savefig(figname, dpi=600)
  plt.close()
  print('Written %s' % figname, file = sys.stdout)

def main():
  outfile = 'result/CDRH3_KD_table_summary_class.csv'
  file_KD = 'result/CDRH3_KD_table_summary.csv'
  file_epi    = 'data/flu_Ab_epi_info.tsv'
  file_seq    = '../Fasta/IGHV1-69_CDRH3.tsv'
  file_flu169 = 'doc/Flu_1-69.xlsx' 
  file_all169 = 'doc/GenBank_1-69.xlsx'
  seq_dict    = seq_to_dict(file_seq)
  epi_dict    = parse_epi_file(file_epi)
  ab_dict_flu = extract_Ab_info(file_flu169)
  ab_dict_all = extract_Ab_info(file_all169)
  ab_name_flu = set([ab_dict_flu[i]['Name'] for i in ab_dict_flu])
  ab_geneID_flu = set([ab_dict_flu[i]['VH Genbank ID'] for i in ab_dict_flu])
  ab_geneID_all = set([ab_dict_all[i]['Genbank ID'] for i in ab_dict_all])
  ab_geneID_other = ab_geneID_all-ab_geneID_flu
  df = process_KD_file(file_KD, ab_name_flu, ab_geneID_flu, ab_geneID_other, epi_dict, seq_dict, outfile)
    
if __name__ == "__main__": 
  main()
