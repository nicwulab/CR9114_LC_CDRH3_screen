# only variable needed to change


PROJECT_PATH='/home/yiquan2/CR9114_CDRH3'
REF = PROJECT_PATH + '/ref/ref.fasta'
merged_FQ = PROJECT_PATH + '/result/{SAMPLENAME}.assembled.fastq'
SAMPLENAMES, = glob_wildcards(merged_FQ)
trimmed_FQ = PROJECT_PATH + '/result/{SAMPLENAME}.trimmed.fastq.gz'
untrimmed_FQ = PROJECT_PATH + '/result/{SAMPLENAME}.untrimmed.fastq.gz'
RESULT_PATH = PROJECT_PATH + '/result/{SAMPLENAME}'
TABLE = RESULT_PATH + '_count.tsv'

rule all:
    input:
        expand(TABLE, SAMPLENAME=SAMPLENAMES),
        expand([trimmed_FQ,untrimmed_FQ], SAMPLENAME=SAMPLENAMES)
rule cutadapt_rm_primer:
    input:
        FQ = merged_FQ
    output:
        trimmed_FQ=trimmed_FQ,
        untrimmed_FQ = untrimmed_FQ
    shell:
        '''cutadapt '''\
        '''-a "^ATCTGAGGACACGGCCGTGTATTACTGT...TGGGGCCAAGGGACCACGGTCACCGTC;e=0.2" -a "^CAGTTTTAGCAGCGGCCCAGCCGGCC...GGACAACCAAAGGCTGCTCCTTCTGT;e=0.2" -a "^CAGTTTTAGCAGCGGCCCAGCCGGCC...AGGACCGTGGCAGCACCTTCCGTGTT;e=0.2" '''\
        '''--untrimmed-output {output.untrimmed_FQ} '''
        '''-O 10 -o {output.trimmed_FQ} {input.FQ}'''

rule fq2count:
    input:
        trimmed_FQ
    params:
        REF_FA=REF
    output:
        TABLE
    shell:
        'python ./fastq2count.py -i {input} -r {params} -o {output}'










