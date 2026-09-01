## Chloe McCreery, Dulguun Amgalan
## 10/23/2023
## Genome Wide Validation Screen
## Demultiplex FASTQ files



########UPDATE
## Directory for gRNA pool validation and alignment: 

sdev -m 32G -t 96:00:00 -p engreitz
conda activate EngreitzLab

PROJECT=/oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/CRISPRiFlowFish/; cd $PROJECT


mkdir -p $PROJECT/{fastq,config,results,log} 

#git clone git@github.com:EngreitzLab/crispri-flowfish.git




################################################################
## Get FASTQ
## Demultiplex

quick-sub -s $PROJECT/log/demultiplex.sh -o $PROJECT/log/demultiplex.out -m 16G -t 08:00:00 \
  "bcl2fastq \
  --runfolder-dir /oak/stanford/groups/engreitz/Projects/SequencingRuns//231020_NB551514_0098_AHTNK5BGXT/ \
  --output-dir $PROJECT/fastq/ \
  --sample-sheet $PROJECT/SampleSheet.csv \
  --no-lane-splitting \
  --create-fastq-for-index-reads \
  --barcode-mismatches 1 \
  --mask-short-adapter-reads 8" 


#Briefly checking fastq contents
#for fastq in fastq/*R1_001.fastq.gz; do echo ${fastq};echo $(zcat ${fastq}|wc -l)/4|bc; done
#for fastq in fastq/*R1_001.fastq.gz; do echo ${fastq}; echo zcat ${fastq} | awk 'NR%4==2' | awk '![x$0]++'; done
#zcat fastq/211119-10-PilotEC-D4-Dox0h-Biorep1-Positive-PCRrep2_S10_R1_001.fastq.gz | awk 'NR%4==2' 


################################################
## Run CRISPRi gRNA quantification pipeline
#git clone git@github.com:EngreitzLab/crispri-flowfish.git

cd fastq  

#Trimmed to 19bp spacers 
for file in *R1_001.fastq.gz; do 
  filename=$(basename "$file")
  output="${filename%.fastq.gz}_trimmed.fastq.gz"
  zcat $filename | cutadapt -u -6 - | gzip > $output
done 

cd ..

BOWTIE_INDEX=$PROJECT/config/2212205_Genome_Wide_Library.fa
bowtie2-build $BOWTIE_INDEX $BOWTIE_INDEX

snakemake \
  -s crispri-flowfish/workflow/Snakefile \
  --configfile config/config.json \
  --cores 1 \
  --jobs 20 \
  --rerun-incomplete \ 
  -k \
  --cluster "sbatch -n 1 -c 1 --mem 8G -t 30:00 -p engreitz -J FF_{rule} -o log/{rule}_{wildcards} -e log/{rule}_{wildcards}"