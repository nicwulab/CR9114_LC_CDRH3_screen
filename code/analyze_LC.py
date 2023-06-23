#!/usr/bin/python
import os
import sys
import logomaker
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

def parse_ab_table(filename):
  datasheet = pd.read_excel(filename, sheet_name='Sheet1')
  Ab_dict = {}
  for index, info_dict in datasheet.iterrows():
    for i in info_dict.keys():
      if 'V_gene' in i:
        info_dict[i] = str(info_dict[i]).rsplit('*')[0]
      elif 'SHM' in i:
        SHMs = [SHM for SHM in str(info_dict[i]).rsplit(',') if SHM[0] not in ['n','s']]
        info_dict[i] = str(','.join(SHMs))
        info_dict[i+'_count'] = str(len(SHMs))
      else:
        info_dict[i] = str(info_dict[i])
    Ab_dict[info_dict['Name']] = info_dict
  return Ab_dict

def count_to_freq(df, colname):
  df[colname+'_freq'] = (df[colname]+1)/(df[colname].sum()+len(df))
  return (df)

def rehashing(Ab_dict, colname):
  new_dict = {}
  for Ab in Ab_dict:
    new_dict[Ab] = Ab_dict[Ab][colname]
  return new_dict

def normalize_score(df, Ab, colname):
  score_non = df.loc[(df['ID'].str.contains('stop')) & (df['input_freq']>0.002)][colname].mean()
  score_WT  = float(df.loc[df['ID']==Ab][colname])
  df[colname] = (df[colname]-score_non)/(score_WT-score_non)
  return df

def process_count_LC(filename, Ab_dict, Abs, epi_dict, chain):
  df  = pd.read_csv(filename)
  df  = df[df['ID'].isin(map(str,Abs.keys()))]
  df  = df.fillna(0)
  df  = df.drop('Unnamed: 0', axis=1)
  for colname in df:
    if colname not in ['ID', 'concentration']:
       df = count_to_freq(df, colname)
  df['input_freq']  = df[chain+'_input_Expression_1_freq']*0.25 + \
                      df[chain+'_input_Expression_2_freq']*0.25 + \
                      df[chain+'_input_Binding_1_freq']*0.25 + \
                      df[chain+'_input_Binding_2_freq']*0.25
  df['exp_score_1'] = np.log10(df[chain+'_bin1_Expression_1_freq']/df[chain+'_bin0_Expression_1_freq'])
  df['exp_score_2'] = np.log10(df[chain+'_bin1_Expression_2_freq']/df[chain+'_bin0_Expression_2_freq'])
  df['bind_score_1'] = np.log10(df[chain+'_bin1_binding_1_freq']/df[chain+'_bin0_binding_1_freq'])
  df['bind_score_2'] = np.log10(df[chain+'_bin1_binding_2_freq']/df[chain+'_bin0_binding_2_freq'])
  for colname in ['exp_score_1', 'exp_score_2', 'bind_score_1', 'bind_score_2']:
    df = normalize_score(df, 'CR9114', colname)
  df['exp_score'] = (df['exp_score_1'] + df['exp_score_2'])/2
  df['bind_score'] = (df['bind_score_1'] + df['bind_score_2'])/2
  df['sequence_aa'] = df['ID'].map(rehashing(Ab_dict, 'VL_AA'))
  df['sequence_nuc'] = df['ID'].map(Abs)
  df['Vgene'] = df['ID'].map(rehashing(Ab_dict, 'Light_V_gene'))
  df['Jgene'] = df['ID'].map(rehashing(Ab_dict, 'Light_J_gene'))
  df['antigen'] = df['ID'].map(rehashing(Ab_dict, 'Antigen'))
  df['SHM_heavy'] = df['ID'].map(rehashing(Ab_dict, 'SHM_heavy'))
  df['SHM_light'] = df['ID'].map(rehashing(Ab_dict, 'SHM_light'))
  df['SHM_heavy_count'] = df['ID'].map(rehashing(Ab_dict, 'SHM_heavy_count'))
  df['SHM_light_count'] = df['ID'].map(rehashing(Ab_dict, 'SHM_light_count'))
  df['CDRL3'] = df['ID'].map(rehashing(Ab_dict, 'CDRL3_AA'))
  df['CDRL3_len'] = df['CDRL3'].str.len()
  df['CDRL3_motif'] = df['CDRL3'].str.slice(start=2,stop=3)+df['CDRL3'].str.slice(start=-2,stop=-1)
  df['epi'] = df['ID'].map(epi_dict) 
  df.to_csv(filename.replace('.csv','_score.csv'),index=False, na_rep='nan')
  return df

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

def write_CDRL3_seq_logo(df, gene, CDRL3_len, high_bind_cutoff, low_bind_cutoff):
  df = df.loc[(df['Vgene']==gene) & (~df['ID'].str.contains('stop'))][['ID','bind_score','exp_score','Jgene','CDRL3']]
  CDRL3_list_high_bind = list(df.loc[df['bind_score']>high_bind_cutoff]['CDRL3'])
  CDRL3_list_low_bind = list(df.loc[df['bind_score']<low_bind_cutoff]['CDRL3'])
  Jgene_list_high_bind = list(df.loc[df['bind_score']>high_bind_cutoff]['Jgene'])
  Jgene_list_low_bind = list(df.loc[df['bind_score']<low_bind_cutoff]['Jgene'])
  make_sequence_logo(CDRL3_list_high_bind, 11, 'graph/seqlogo_'+gene+'_high_bind.png')
  make_sequence_logo(CDRL3_list_low_bind, 11, 'graph/seqlogo_'+gene+'_low_bind.png')

def write_CDRL3_motif_logo(df, high_bind_cutoff, low_bind_cutoff):
  motif_list_high_bind = list(df.loc[df['bind_score']>high_bind_cutoff]['CDRL3_motif'])
  motif_list_low_bind = list(df.loc[df['bind_score']<low_bind_cutoff]['CDRL3_motif'])
  print ("# of variant with binding score > %f: %i" % (high_bind_cutoff, len(motif_list_high_bind)))
  print ("# of variant with binding score < %f: %i" % (low_bind_cutoff, len(motif_list_low_bind)))
  make_sequence_logo(motif_list_high_bind, 2, 'graph/seqlogo_motif_high_bind.png')
  make_sequence_logo(motif_list_low_bind, 2, 'graph/seqlogo_motif_low_bind.png')

def analyze_SHM(df, bind_cutoff):
  df = df.loc[(~df['ID'].str.contains('stop')) & \
              (df['input_freq']>0.002) & \
              (df['Vgene'].str.contains('IGLV1')) & \
              ((~df['CDRL3_motif'].str.slice(start=0,stop=1).isin(['F','Y','W'])) | \
              (df['CDRL3_motif'].str.slice(start=1,stop=2).isin(['F','Y','W'])))]
  SHM_high_bind = set((','.join(list(df.loc[df['bind_score']>bind_cutoff]['SHM_light']))).rsplit(','))
  SHM_low_bind  = set((','.join(list(df.loc[df['bind_score']<bind_cutoff]['SHM_light']))).rsplit(','))
  print (df.loc[df['bind_score']>bind_cutoff])
  print (SHM_high_bind - SHM_low_bind)

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

def main():
  outfile      = 'result/LC_count_table_epi.csv'
  LC_file      = 'result/LC_count_table.csv'
  file_ab      = 'doc/paired 1-69 antibodies_shm.xlsx'
  file_epi     = 'data/flu_Ab_epi_info.tsv'
  LC_ref_file  = 'Fasta/lib_LC.fa'
  LC_Abs      = SeqIO.to_dict(SeqIO.parse(LC_ref_file, "fasta"))
  for ID in LC_Abs.keys():
    LC_Abs[ID] = str(LC_Abs[ID].seq)
  Ab_dict     = parse_ab_table(file_ab)
  epi_dict    = parse_epi_file(file_epi)
  df = process_count_LC(LC_file, Ab_dict, LC_Abs, epi_dict, 'L')
  df = df.loc[(~df['ID'].str.contains('stop')) & \
              (df['input_freq']>0.002)]
  write_CDRL3_motif_logo(df, 0.8, 0.8)

if __name__ == "__main__":
  main()
