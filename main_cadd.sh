#!/usr/bin/env bash

sbatch_fname="main_array_cadd.sbatch"
project="vep"
lista="lista_archivos.txt" # Asegúrate de que este nombre coincide con el que tienes
mkdir -p /home/ntoro/clinvar/clinvar_cortado/jobs_vep


# Contar líneas (tienes unos 30 archivos en esa lista)
n_files=$(wc -l < $lista)

# Lanzar array (limite de 10 a la vez por ejemplo)
sbatch -a 1-${n_files} $sbatch_fname
