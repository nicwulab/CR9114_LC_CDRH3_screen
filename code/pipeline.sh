#!/bin/bash

snakemake -s ./pipeline_merge.smk -j 64
snakemake -s ./pipeline_count.smk -j 64


