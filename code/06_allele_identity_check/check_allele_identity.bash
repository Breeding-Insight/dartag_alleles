#!/bin/bash
#SBATCH --job-name=Mummer_coordi      # Job name
#SBATCH --output=Mummer_coordi_%j.out  # Output log file (with job ID)
#SBATCH --error=Mummer_coordi_%j.err   # Error log file (with job ID)
#SBATCH --partition=regular                   # Partition name
#SBATCH --nodes=1                             # Number of nodes
#SBATCH --ntasks=1                            # Number of tasks (typically 1 per script)
#SBATCH --cpus-per-task=1                     # Number of CPUs per task
#SBATCH --mem=500G                            # Total memory for the job
#SBATCH --mail-user=dz359@cornell.edu        # Your email for notifications
#SBATCH --mail-type=ALL                       # Notifications for job BEGIN, END, FAIL
#SBATCH --time=10-00:00:00                      # Maximum time for the job (10 days)

# map alleles to reference genomes
bwa mem \
    -t 1 \
    -M \
    -R "@RG\tID:cranberry\tSM:cranberry\tPL:illumina" \
    /path/to/reference_genome.fasta \
    /path/to/allele_sequences.fastq.gz | \
