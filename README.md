## Probing the sequence constraints of CR9114 CDR H3 and light chain

## ENV setup
```
conda env create -f environment.yml
```

or create ENV with commands

```
conda create -n CDRH3 python=3.9
conda activate CDRH3
conda install -c bioconda -c anaconda -c conda-forge \
  pear \
  biopython \
  cutadapt \
  snakemake \
  argparse \
  logomaker \
  openpyxl 
```

## Pipeline Usage

Run the snakemake pipeline to get the sequence read count:
```
./pipeline.sh
```

for KD NLS regression, run the notebook [CDRH3_Tite_seq_plot.Rmd](./code/CDRH3_Tite_seq_plot.Rmd) by R

## Analysis of Light chain variants
1. Run ``python3 code/analyze_LC.py`` to compute expression and binding scores for light chain variants and plot sequence logos.
2. Run ``Rscript code/plot_LC.R`` to plot analysis of light chain variants.

## Analysis of CDR H3 variants
