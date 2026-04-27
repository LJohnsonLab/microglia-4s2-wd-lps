#!/bin/bash 
#SBATCH --job-name=nextflow_rnaseq
#SBATCH --output=nextflow_%j.out 
#SBATCH --error=nextflow_%j.err 
#SBATCH --time=24:00:00 
#SBATCH --cpus-per-task=1 
#SBATCH --mem=4G
#SBATCH --account=coa_lajo247_uksr
#SBATCH --partition=normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user jm.arbones@uky.edu

# Load nextflow
module load ccs/conda/nextflow/24.10.4
# Load the environment  
export SCRATCH=/scratch/jar301
export NXF_SINGULARITY_CACHEDIR=$SCRATCH/.singularity/cache
export SINGULARITY_CACHEDIR=$SCRATCH/.singularity/cache

# Run the pipeline in nextflow
nextflow run nf-core/rnaseq \
--input /scratch/jar301/TR/SRP453154/20250731-sample_sheet_SRP453154.csv \
--outdir /scratch/jar301/TR/SRP453154/rnaseq_results \
--aligner star_salmon \
--fasta /home/jar301/reference_genome/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz \
--gtf /home/jar301/reference_genome/Mus_musculus.GRCm39.114.gtf.gz \
--skip_markduplicates \
--igenomes_ignore \
--genome null \
--email jm.arbones@uky.edu \
-profile singularity \
-c /scratch/jar301/nextflow.config
