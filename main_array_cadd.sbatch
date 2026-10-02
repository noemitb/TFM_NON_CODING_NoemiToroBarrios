#!/usr/bin/env bash

#SBATCH -c 4         
#SBATCH --mem=8G
#SBATCH --output=jobs_vep/%A_%a.out
#SBATCH --error=jobs_vep/%A_%a.err

# Noemí Toro Barrios | TFM - Revisión de herramientas y recursos para la anotación y priorización de variantes en regiones no codificantes
#Anotación de variantes con CADD mediante VEP

#Cargo el módulo de singularity 

module load singularity

#Directorio donde residen los VCFs de entrada y nombre de archivos.
VCF_DIR="/home/ntoro/clinvar/clinvar_cortado/split"
LISTA="lista_archivos.txt"

line=$SLURM_ARRAY_TASK_ID

#Obtenego el ID
NAME_ID=$(sed -n "${line}p" ${LISTA})

echo "Procesando ID: $NAME_ID"


#Defino solo el nombre del archivo (sin ruta completa) para usarlo dentro del contenedor
INPUT_FILENAME="${NAME_ID}.vcf.gz"
OUTPUT_FILENAME="${NAME_ID}_anotado.vcf.gz" # Añadido .gz porque usas bgzip


#ESTRUCTURA DE MONTAJES (-B):
#Monto el directorio de entrada (VCF_DIR) en /input_dir
#Monto el directorio actual ($PWD) en /output_dir para guardar resultados
#Monto caches y plugins

singularity exec \
    -B "$VCF_DIR":/input_dir \
    -B "$PWD"/cadd:/output_dir \
    -B /mnt/lustre/scratch/CBRA/data/vep_cache/:/cache_dir \
    -B /mnt/lustre/home/ntoro/:/fasta_dir \
    -B /mnt/lustre/scratch/CBRA/data/vep_cache/Plugins/resources/CADD_v1_7_grch38:/plugin_dir \
    docker://ensemblorg/ensembl-vep \
    vep \
    --input_file /input_dir/"$INPUT_FILENAME" \
    --output_file /output_dir/"$OUTPUT_FILENAME" \
    --vcf \
    --compress_output bgzip \
    --force_overwrite \
    --fork 4 \
    --assembly GRCh38 \
    --cache --offline \
    --dir_cache /cache_dir \
    --fasta /fasta_dir/hg38.fa \
    --custom /whole_genome_SNVs.tsv.gz,/plugin_dir/gnomad.genomes.r4.0.indel.tsv.gz
