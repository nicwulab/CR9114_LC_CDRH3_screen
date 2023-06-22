## Probing the sequence constraints of CR9114 CDR H3 and light chain

## ENV setup
```
conda env create -f env.yml
```

or create ENV with commands

```conda create -n CDRH3 python=3.9
conda activate CDRH3
conda install -c bioconda -c anaconda -c conda-forge \
  pear \
  biopython \
  cutadapt \
  snakemake \
  argparse
```

## Pipeline Usage

Run the snakemake pipeline to get the sequence read count:
```
./pipeline.sh
```

for KD NLS regression, run the notebook [CDRH3_Tite_seq_plot.Rmd](./code/CDRH3_Tite_seq_plot.Rmd) by R
