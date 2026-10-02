#!/usr/bin/env python3

# Noemí Toro Barrios | TFM: Revisión de herramientas y recursos para la anotación y priorización de variantes en regiones no codificantes.

# Script para extraer la anotación de PDIVAS, clasificar las variantes usando el umbral 0.082 y añadir el resultado a la tabla de splicing.

# Cargo todas las librerías, incluido el archivo donde tengo definidas las funciones pero sin el .py

import sys
import gzip
import os
import utils


# Tabla que voy a utilizar como referencia y nombre de la tabla final

tabla_entrada = "megatabla_resultados_ADA.tsv"
tabla_salida = "megatabla_resultados_splicing.tsv"


# Parseo la salida de PDIVAS para obtener el score de cada variante.

def obtener_scores_pdivas(pdivas_str):  # La anotación de la herramienta aparece en el campo PDIVAS.

    # Compruebo que exista anotación antes de intentar sacar el score.

    if not pdivas_str or pdivas_str == ".":

        return -1.0  # Utilizo -1 para indicar después que la variante no ha sido anotada.

    else:

        anotaciones = pdivas_str.split(",")

        scores = []

        # Puede haber varias anotaciones para la misma variante, así que las recorro todas.

        for anot in anotaciones:

            partes = anot.split("|")

            # Si no tengo suficientes campos no puedo recuperar el score y salto esta anotación.

            if len(partes) <= 1:

                continue

            score_raw = partes[1]

            # Intento convertir el score a número. Si no puedo, simplemente salto ese valor.

            try:

                scores.append(float(score_raw))

            except ValueError:

                continue

        # Si tengo más de un score, me quedo con el valor máximo.

        if scores:

            return max(scores)

        else:

            return -1.0


def main():

    # Compruebo que he indicado el archivo con la lista de VCFs que quiero procesar.

    if len(sys.argv) < 2:

        print("No has indicado lista de VCFs para extraer la anotación")
        sys.exit(1)

    # Cargo el archivo donde tengo las rutas de los VCFs.

    archivo_lista = sys.argv[1]

    lista_rutas = []

    with open(archivo_lista, "r") as f_lista:  # Abro el archivo en modo lectura.

        for linea in f_lista:

            ruta = (
                linea.strip()
            )  # Quito posibles espacios o saltos de línea.

            if ruta:

                lista_rutas.append(
                    ruta
                )  # Si la línea contiene una ruta, la añado a mi lista.

    # Creo un diccionario donde voy a guardar el score de PDIVAS de cada variante.

    memoria_pdivas = {}

    for indice, ruta_archivo in enumerate(lista_rutas, start=1):

            nombre_archivo = os.path.basename(
                ruta_archivo
            )  # Me quedo solamente con el nombre del archivo.

            print(f"Leyendo y extrayendo datos de: {nombre_archivo}")

            if not os.path.exists(ruta_archivo):  # Compruebo que el archivo existe.

                print(f"No se ha encontrado el archivo {ruta_archivo}")
                continue

            # Guardo los identificadores de las variantes que ya aparecen en mi tabla de entrada.

            variantes_spliceai = set()

            # Me quedo con las variantes de la tabla ya generada porque son las que han pasado los filtros anteriores.

            with open(tabla_entrada, "r", encoding="utf-8") as f_in:

                # La primera línea es la cabecera, así que la salto.

                next(f_in)

                for linea in f_in:

                    variante_id = linea.split("\t")[0]

                    variantes_spliceai.add(variante_id)

                # Abro el VCF anotado con PDIVAS y recorro todas sus variantes.

                with open(ruta_archivo, "r", encoding="utf-8") as f_vcf:

                    for variant in f_vcf:

                        if variant.startswith(
                        "#"
                        ):  # Las líneas que empiezan por # pertenecen a la cabecera y no las necesito.

                            continue

                        if not variant:  # Si la línea está vacía, la salto.

                            continue

                        cols = variant.split("\t")

                        # Quito "chr" para mantener el mismo formato de identificador que utilizo en mis tablas.

                        cromosoma = cols[0].replace("chr", "")

                        # Formo el identificador de la variante para poder relacionarla con la tabla de entrada.

                        variante_id = f"{cromosoma}:{cols[1]}:{cols[3]}:{cols[4]}"

                        datos_info = utils.parse_info(cols[7])

                        pdivas_anot = datos_info.get("PDIVAS", ".")

                        # Obtengo el score de PDIVAS y lo guardo asociado a la variante.

                        pdivas_score = obtener_scores_pdivas(pdivas_anot)

                        memoria_pdivas[variante_id] = pdivas_score

    # Abro la tabla de entrada y creo la nueva tabla donde añadiré la clasificación de PDIVAS.

    with open(tabla_entrada, "r", encoding="utf-8") as f_in, open(tabla_salida, "w", encoding="utf-8") as f_out:

        cabecera_original = f_in.readline().strip()  # Leo la cabecera y le quito el salto de línea.

        # Añado la nueva columna de PDIVAS a la cabecera.

        f_out.write(f"{cabecera_original}\tPDIVAS(0.082)\n")

        for linea in f_in:

            linea_limpia = linea.strip()  # Quito el salto de línea.

            columnas_tsv = linea_limpia.split("\t")  # Separo las columnas por tabuladores.

            variante_id = columnas_tsv[0]  # La primera columna contiene la variante. Ejemplo: 14:55466853:A:G

            # Busco la variante en el diccionario donde he guardado los scores.
            # Si no la encuentro, utilizo -1 para marcarla como No anotada.

            scores = memoria_pdivas.get(variante_id, -1.0)

            # Clasifico la variante como patogénica o benigna utilizando el umbral 0.082.

            if scores == -1:

                class_pdivas = "No anotada"

            elif scores >= 0.082:

                class_pdivas = "Patogénica"

            else:

                class_pdivas = "Benigna"

            # Escribo la línea original y añado al final la clasificación de PDIVAS.

            f_out.write(f"{linea_limpia}\t{class_pdivas}\n")


if __name__ == "__main__":
    main()