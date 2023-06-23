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
  openpyxl 
```

## Pipeline Usage

Run the snakemake pipeline to get the sequence read count:
```
./pipeline.sh
```

for KD NLS regression, run the notebook [CDRH3_Tite_seq_plot.Rmd](./code/CDRH3_Tite_seq_plot.Rmd) by R

## Analysis of Light chain variants
1. Run ``python3 code/analyze_LC.py`` to compute expression and binding scores for light chain variants and plot sequence logos using [Logomaker](https://logomaker.readthedocs.io/en/latest/).
2. Run ``Rscript code/plot_LC.R`` to plot analysis of light chain variants.

## Analysis of CDR H3 variants
1. Run ``python3 code/process_KD_table.py`` to process the KD results.
2. Run ``python3 code/compute_expression_score.py`` to compute expression scores for CDR H3 variants.
3. Run ``Rscript code/plot_CDRH3_QC.R`` to plot the correlation between replicates.
4. Run ``Rscript code/plot_CDRH3_KD_distrib.R`` to plot the distribution of expression and binding scores, and analyzed the effect of Y98 on binding.
5. Run ``Rscript code/plot_CDRH3_exp_vs_KD.R`` to plot expression scores vs binding scores.
6. Run ``Rscript code/plot_CDRH3_KD_heatmap.R`` to plot the heatmap to show the effect of single mutants on binding.
