#!/usr/bin/python
import sys
from Bio import SeqIO
from collections import Counter
import os
import argparse
import gzip

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
    "---":"-"}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i + 3].upper()
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def sum_mut(aa1,aa2):
    return sum ( aa1[i] != aa2[i] for i in range(len(aa1)) )


def call_mutid(mutpep,refseq,shift):
  mut_id_ls = []
  assert (len(mutpep) == len(refseq))
  for n in range(len(mutpep)):
    pos = n+shift
    if refseq[n]!=mutpep[n]:
       mut_id_ls.append(refseq[n]+str(pos)+mutpep[n])
  return mut_id_ls



def cal_ref_dict(ref):
    # create a reference dictionary
    Rrecords = SeqIO.parse(ref, "fasta")
    Rs=[]
    ref_dict={}
    for record in Rrecords:
        seq = str(record.seq)
        id = str(record.id)
        if id in ref_dict.values():
            record.id = record.id + '_1'
            id = id + '_1'
        ref_dict[seq]=id
        Rs.append(record)
    SeqIO.write(Rs, "/ref/ref_clean.fasta", "fasta")
    return ref_dict

def cal_fastq_dic(fastq,ref_dict):

    print ("reading %s" % fastq)
    with gzip.open(fastq,'rt') as handle:
        Rrecords = SeqIO.parse(handle,"fastq")

        id_count_ls = []
        mut_ls = []
        error_read=0
        for record in Rrecords:
            seq = str(record.seq)
            #filter out ambiguous read
            if 'N' in seq:continue
            if seq in ref_dict.keys():
                id = ref_dict[seq]
                id_count_ls.append(id)
            else:
                mut_ls.append(seq)
                error_read+=1

        print('incorrect barcode, a total of %s reads were excluded' %error_read)
        lib_dict = Counter(id_count_ls)
        mut_dict = Counter(mut_ls)
    return lib_dict,mut_dict

def write_mut_table(mut_dic,outfilename):
  outfile = open(outfilename,'w')
  outfile.write("\t".join(['Mutation', 'count'])+"\n")
  for mut in mut_dic.keys():
    Mutation = mut
    count = mut_dic[mut]
    outfile.write("\t".join(map(str,[Mutation, count]))+"\n")
  outfile.close()
  print('Written %s' %outfilename)

def write_count_table(id_count_dic,outfilename):
  outfile = open(outfilename,'w')
  outfile.write("\t".join(['ID', 'count'])+"\n")
  for id in id_count_dic.keys():
    count = id_count_dic[id]
    outfile.write("\t".join(map(str,[id, count]))+"\n")
  outfile.close()
  print('Written %s' %outfilename)

def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-r', '--reference_fasta', type=str,  help="reference fasta file path.")
    parser.add_argument('-i', '--input', type=str,  help="fastq file path.")
    parser.add_argument('-o', '--output', type=str,  help="output count tsv file path.")


    args = parser.parse_args()

    ref_dict = cal_ref_dict(args.reference_fasta)
    lib_dict,mut_dict = cal_fastq_dic(args.input,ref_dict)
    write_count_table(lib_dict, args.output)
    mut_output=args.output.split('.')[0] + '_mut.tsv'
    write_mut_table(mut_dict, mut_output)


if __name__ == "__main__":
  main()

