#!/bin/bash
#
#SBATCH --job-name=quick-sub
#SBATCH --ntasks=1
#SBATCH --partition=engreitz,normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=2:00:00
#SBATCH --output /oak/stanford/groups/engreitz/Users/chloemc/231023_GenomeWide//log/demultiplex.out
bcl2fastq --runfolder-dir /oak/stanford/groups/engreitz/Projects/SequencingRuns//231020_NB551514_0098_AHTNK5BGXT/ --output-dir /oak/stanford/groups/engreitz/Users/chloemc/231023_GenomeWide//fastq/ --sample-sheet /oak/stanford/groups/engreitz/Users/chloemc/231023_GenomeWide//SampleSheet_swapped.csv --no-lane-splitting --create-fastq-for-index-reads --barcode-mismatches 1 --mask-short-adapter-reads 8
