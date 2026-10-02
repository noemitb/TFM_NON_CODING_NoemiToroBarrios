#!/usr/bin/env bash

#SBATCH -c 4
#SBATCH --mem=128G
#SBATCH --output=jobs/%A_%a.out
#SBATCH --error=jobs/%A_%a.err

# Noemí Toro Barrios | TFM - Revisión de herramientas y recursos para la anotación y priorización de variantes en regiones no codificantes
#Anotación de variantes con PDIVAS

module load anaconda3

source $CONDA_PREFIX/etc/profile.d/conda.sh

export PATH="/mnt/lustre/home/ntoro/.conda/envs/PDIVAS_FINAL/bin:$PATH"

line=$SLURM_ARRAY_TASK_ID

module load singularity
module load bcftools/1.9

#Defino las rutas y archivos que voy a utilizar

VCF_DIR="/home/ntoro/clinvar/clinvar_cortado/split/troceado"
LISTA="lista_archivos.txt"
CHR_MAP="chr_map.txt"

#Obtengo el ID

NAME_ID=$(sed -n "${line}p" ${LISTA})

echo "Procesando ID: $NAME_ID"

#Defino los nombres de los archivos para usarlos dentro del contenedor

INPUT_FILENAME_ORIGINAL="${NAME_ID}.vcf.gz"
INPUT_FILENAME_PROCESADO="${NAME_ID}_procesado.vcf.gz"
OUTPUT_FILENAME_VEP="${NAME_ID}_anotado_vep.vcf.gz"

#Modifico el chr en los vcfs

bcftools annotate --rename-chrs "$CHR_MAP" "$VCF_DIR/$INPUT_FILENAME_ORIGINAL" -O u | \
bcftools norm -m -any -O u | \
bcftools filter -i 'STRLEN(REF)<50 && STRLEN(ALT)<50' -O u | \
bcftools view -e 'ALT="*"' -O z -o "$VCF_DIR/$INPUT_FILENAME_PROCESADO"

tabix -p vcf "$VCF_DIR/$INPUT_FILENAME_PROCESADO"

#Monto las rutas necesarias dentro del contenedor

singularity exec \
    -B "$VCF_DIR":/input_dir \
    -B "$PWD":/output_dir \
    -B /mnt/lustre/scratch/CBRA/data/vep_cache/:/cache_dir \
    -B /mnt/lustre/home/ntoro/:/fasta_dir \
    -B /home/ntoro/clinvar/clinvar_cortado/vep_files/:/plugin_dir \
    docker://ensemblorg/ensembl-vep \
    vep \
    --input_file /input_dir/"$INPUT_FILENAME_PROCESADO" \
    --output_file /output_dir/"$OUTPUT_FILENAME_VEP" \
    --vcf \
    --compress_output bgzip \
    --force_overwrite \
    --fork 4 \
    --assembly GRCh38 \
    --cache --offline \
    --cache_version 115 \
    --dir_cache /cache_dir \
    --fasta /fasta_dir/hg38.fa \
    --hgvs \
    --pick_allele_gene \
    --custom /plugin_dir/ConSplice.50bp_region.inverse_proportion_refo_hg38.bed.gz,ConSplice,bed,overlap,0 \
    --plugin MaxEntScan,/plugin_dir/fordownload,SWA,NCSS \
    --fields "Consequence,SYMBOL,Gene,INTRON,HGVSc,STRAND,ConSplice,MES-SWA_acceptor_diff,MES-SWA_acceptor_alt,MES-SWA_donor_diff,MES-SWA_donor_alt"

conda deactivate

tabix -p vcf "$OUTPUT_FILENAME_VEP"

export PATH="/mnt/lustre/home/ntoro/.conda/envs/PDIVAS_FINAL/bin:$PATH"

#Defino SpliceAI y la referencia

SPLICEAI="/mnt/lustre/home/ntoro/.conda/envs/PDIVAS_FINAL/bin/spliceai"
REF_DIR="/mnt/lustre/home/ntoro/hg38.fa"
REF="grch38"

#Defino los archivos de entrada y salida de SpliceAI

INPUT_FILE_SPLICEAI="./${OUTPUT_FILENAME_VEP}"
OUTPUT_FILE_SPLICEAI="${NAME_ID}_anotado_spliceai.vcf"

echo "Anotando archivo: $INPUT_FILE_SPLICEAI"

#Ejecuto SpliceAI con los parámetros necesarios para PDIVAS

$SPLICEAI \
    -I "$INPUT_FILE_SPLICEAI" \
    -O "$OUTPUT_FILE_SPLICEAI" \
    -R "$REF_DIR" \
    -A "$REF" \
    -D 300 \
    -M 1

echo "Guardado en: $OUTPUT_FILE_SPLICEAI"

conda deactivate

#Defino el archivo de salida de PDIVAS

OUTPUT_FILE_PDIVAS="${NAME_ID}_anotado_pdivas.vcf"

#Ejecuto PDIVAS

pdivas predict \
    -I "$OUTPUT_FILE_SPLICEAI" \
    -O "$OUTPUT_FILE_PDIVAS" \
    -F off
