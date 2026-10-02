#!/usr/bin/env bash

sbatch_fname="main_array_pdivas.sbatch"
project="pdivas"
lista="lista_archivos.txt" # Asegúrate de que este nombre coincide con el que tienes
mkdir -p /home/ntoro/clinvar/clinvar_cortado/pdivas/jobs


# Contar líneas (tienes unos 30 archivos en esa lista)
n_files=$(wc -l < $lista)

# Lanzar array (limite de 10 a la vez por ejemplo)
sbatch -a 1-${n_files}%30 $sbatch_fname
