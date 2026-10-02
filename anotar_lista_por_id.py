#!/usr/bin/env python3

# Noemí Toro Barrios | TFM: Revisión de herramientas y recursos para la anotación y priorización de variantes en regiones no codificantes.

# Script para añadir a los VCFs la puntuación LOEUF de gnomAD utilizando el símbolo del gen o su identificador Ensembl.

import os
import gzip
import re


# Rutas de los archivos que voy a utilizar

archivo_gnomad = "/mnt/lustre/home/ntoro/clinvar/clinvar_cortado/vep_files/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz"
lista_archivos = "/mnt/lustre/home/ntoro/clinvar/clinvar_cortado/clinvar_vep/lista_vcfs.txt"
dir_salida = "/mnt/lustre/home/ntoro/clinvar/clinvar_cortado/clinvar_vep/anotados_loeuf"


# Creo un diccionario donde voy a guardar el valor LOEUF asociado tanto al símbolo del gen como a su identificador Ensembl.

loeuf_dict = {}

with gzip.open(archivo_gnomad, "rt") as f:

    for line in f:

        # Salto la cabecera del archivo de gnomAD.

        if line.startswith("#gene") or line.startswith("gene\t"):

            continue

        partes = line.strip().split('\t')

        # Compruebo que la línea tenga todas las columnas necesarias y que el valor LOEUF esté disponible.

        if len(partes) > 63 and partes[29] != "NA":

            symbol = partes[0]

            score = partes[29]

            gene_id = partes[63]

            # Guardo el mismo score usando como clave tanto el símbolo del gen como el ID de Ensembl.

            loeuf_dict[gene_id] = score

            loeuf_dict[symbol] = score


# Compruebo que exista el archivo con la lista de VCFs antes de continuar.

if not os.path.exists(lista_archivos):

    exit(1)


# Cargo todas las rutas de los VCFs que quiero anotar.

with open(lista_archivos, 'r') as f:

    rutas = f.read().splitlines()


# Recorro cada uno de los VCFs de la lista.

for vcf_entrada in rutas:

    vcf_entrada = vcf_entrada.strip()

    if not vcf_entrada:

        continue

    nombre_base = os.path.basename(vcf_entrada).replace('.vcf.gz', '')

    # Creo el nombre del VCF de salida manteniendo el nombre del archivo original.

    vcf_final = os.path.join(
        dir_salida,
        f"{nombre_base}_LOEUF_FINAL.vcf.gz"
    )

    # Si el archivo no existe, lo salto y continúo con el siguiente.

    if not os.path.exists(vcf_entrada):

        continue

    with gzip.open(vcf_entrada, "rt") as fin, gzip.open(vcf_final, "wt") as fout:

        for line in fin:

            # Mantengo las líneas de cabecera del VCF y añado la definición del nuevo campo LOEUF_score.

            if line.startswith("#"):

                if line.startswith("#CHROM"):

                    fout.write(
                        '##INFO=<ID=LOEUF_score,Number=1,Type=String,Description="gnomAD v2.1.1 LOEUF score (oe_lof_upper) matched by Gene Symbol or Ensembl ID">\n'
                    )

                fout.write(line)

            else:

                # Quito el salto de línea para poder trabajar con las columnas del VCF.

                line_clean = line.rstrip('\n')

                columnas = line_clean.split('\t')

                # Busco todos los identificadores Ensembl que puedan aparecer en la variante.

                genes_a_comprobar = set(
                    re.findall(r'ENSG\d+', line_clean)
                )

                # Compruebo también si existe el símbolo del gen dentro del campo GENEINFO.

                m_geneinfo = re.search(
                    r'GENEINFO=([^:;]+)',
                    line_clean
                )

                if m_geneinfo:

                    genes_a_comprobar.add(
                        m_geneinfo.group(1)
                    )

                # Busco también los símbolos de genes que aparecen dentro de la anotación CSQ.

                m_csq = re.findall(
                    r'\|([A-Za-z0-9orf\-\.]+)\|ENSG\d+',
                    line_clean
                )

                for sym in m_csq:

                    if sym:

                        genes_a_comprobar.add(sym)

                # Recorro todos los genes encontrados hasta localizar uno que tenga valor LOEUF.

                score_final = None

                for g in genes_a_comprobar:

                    if g in loeuf_dict:

                        score_final = loeuf_dict[g]

                        break

                # Si encuentro un valor LOEUF, lo añado a la columna INFO del VCF.

                if score_final:

                    columnas[7] = (
                        columnas[7]
                        + f";LOEUF_score={score_final}"
                    )

                # Vuelvo a montar la línea del VCF y la escribo en el archivo de salida.

                fout.write(
                    '\t'.join(columnas)
                    + '\n'
                )
