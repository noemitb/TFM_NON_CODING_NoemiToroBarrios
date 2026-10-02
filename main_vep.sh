#!/usr/bin/env bash

sbatch_fname="main_array_vep.sbatch"
lista="lista_archivos.txt"

# Carpeta para los logs de VEP
mkdir -p /home/ntoro/clinvar/clinvar_cortado/jobs_vep

# Contar el número de archivos de la lista
n_files=$(wc -l < "$lista")

echo "Número de archivos: $n_files"
echo "Lanzando: $sbatch_fname"

# Lanzar array, máximo 10 tareas simultáneas
sbatch -a 1-${n_files}%10 "$sbatch_fname"
