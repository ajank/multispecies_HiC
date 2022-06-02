#!/bin/bash

# Run the pipeline locally, using 32 cores
snakemake --snakefile src/snakemake/evoHiC.Snakefile --use-conda -j 32 -k
