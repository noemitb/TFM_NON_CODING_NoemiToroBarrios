#!/usr/bin/env bash

# Noemí Toro Barrios | TFM: Revisión de herramientas y recursos para la anotación y priorización de variantes en regiones no codificantes.

# Script principal para lanzar mediante un array de SLURM la anotación de los VCFs de ClinVar con SpliceAI.

# Indico el script que contiene la ejecución de SpliceAI.

sbatch_fname="main_array_spliceai.sbatch"

project="spliceai"

# Archivo donde tengo los nombres de todos los VCFs que quiero procesar.

lista="lista_archivos.txt"

# Creo la carpeta donde se guardarán los archivos de salida y error de SLURM.

mkdir -p jobs


# Cuento cuántos VCFs tengo en la lista para saber cuántas tareas tiene que tener el array.

n_files=$(wc -l < "$lista")


# Lanzo una tarea del array por cada archivo de la lista.

sbatch -a 1-${n_files} "$sbatch_fname"
