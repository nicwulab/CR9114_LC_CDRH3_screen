# only variable needed to change


PROJECT_PATH='/home/yiquan2/CR9114_CDRH3'

R1 = PROJECT_PATH + '/fastq/{SAMPLENAME}_L001_R1_001.fastq.gz'
R2 = PROJECT_PATH + '/fastq/{SAMPLENAME}_L001_R2_001.fastq.gz'
SAMPLENAMES, = glob_wildcards(R1)

merged_FQ = PROJECT_PATH + '/result/{SAMPLENAME}'



rule all:
    input:
        expand(merged_FQ, SAMPLENAME=SAMPLENAMES),




rule pear_merge_reads:
    input:
        FQ1 = R1,
        FQ2 = R2
    output:
        merged_FQ
    shell:
        'pear -f {input.FQ1} -r {input.FQ2} -o {output} -j 64'


