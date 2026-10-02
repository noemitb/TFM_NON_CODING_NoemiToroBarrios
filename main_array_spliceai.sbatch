#!/usr/bin/env bash

#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH --output=jobs/%A_%a.out
#SBATCH --error=jobs/%A_%a.err

#Script para lanzar un array para la anotación de los VCFs divididos de Clinvar. 
#TFM_NOEMITORO


#Establezco variantes
VCF_DIR="/home/ntoro/clinvar/clinvar_cortado/split"
SPLICEAI="/home/ntoro/.conda/envs/spliceai_env/bin/spliceai" #Para no poner la ruta en la ejecución de la herramienta
LISTA="lista_archivos.txt" #Contienen los nombres para que llame a los VCFs utilizados
REF_DIR="/mnt/lustre/home/ntoro/hg38.fa"
REF="grch38" 



line=$SLURM_ARRAY_TASK_ID

#Leo la lista
NAME_ID=$(sed -n "${line}p" ${LISTA})

echo "Procesando ID: $NAME_ID"

#Establezco los archivos de entradas y salidas la ejecución del programa
INPUT_FILE="${VCF_DIR}/${NAME_ID}.vcf.gz"

OUTPUT_FILE="${NAME_ID}_anotado.vcf"


echo "Anotando archivo: $INPUT_FILE"

#Ejecuto la herramienta con las variables anteriores

$SPLICEAI -I "$INPUT_FILE" -O "$OUTPUT_FILE" -R "$REF_DIR" -A "$REF"


echo "Guardado en: $OUTPUT_FILE"
