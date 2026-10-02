#!/usr/bin/env python3

# Noemí Toro Barrios | TFM: Revisión de herramientas y recursos para la anotación y priorización de variantes en regiones no codificantes.

# Script para extraer las puntuaciones ADA y RF de dbscSNV, clasificarlas usando el umbral 0.6 y añadirlas a la tabla de resultados.

# Cargo todas las librerías, incluido el archivo donde tengo definidas las funciones pero sin el .py

import sys
import gzip
import os
import utils


# Tabla que voy a utilizar como entrada y nombre de la tabla final.

tabla_entrada = "megatabla_resultados_herramientas.tsv"
tabla_salida = "megatabla_resultados_ADA.tsv"


# Parseo la salida de VEP para obtener los scores ADA y RF de dbscSNV.

def obtener_scores_dbscsnv(csq_str):  # La anotación de la herramienta aparece en el campo CSQ.

    # Compruebo que exista anotación antes de intentar sacar los scores.

    if not csq_str or csq_str == ".":

        return -1.0, -1.0  # Devuelvo -1 para ADA y RF si la variante no está anotada.

    adas = []

    rfs = []

    # Como puede haber diferentes anotaciones por transcrito, separo cada uno de ellos.

    transcritos = csq_str.split(',')

    for t in transcritos:

        partes = t.split('|')

        # Si no tengo suficientes campos no puedo recuperar los scores y salto este transcrito.

        if len(partes) <= 26:

            continue

        ada = partes[25]

        rf = partes[26]

        # Si existe score ADA, intento convertirlo a número y lo guardo.

        if ada and ada != ".":

            try:

                adas.append(float(ada))

            except ValueError:

                pass

        # Hago lo mismo con el score RF.

        if rf and rf != ".":

            try:

                rfs.append(float(rf))

            except ValueError:

                pass

    # Si hay varios scores por transcrito, me quedo con el valor máximo de cada predictor.

    max_ada = max(adas) if adas else -1.0

    max_rf = max(rfs) if rfs else -1.0

    return max_ada, max_rf


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

    # Creo un diccionario donde voy a guardar los scores ADA y RF de cada variante.

    memoria_dbscsnv = {}

    for indice, ruta_archivo in enumerate(lista_rutas, start=1):

            nombre_archivo = os.path.basename(
                ruta_archivo
            )  # Me quedo solamente con el nombre del archivo.

            print(f"Leyendo y extrayendo datos de: {nombre_archivo}")

            if not os.path.exists(ruta_archivo):  # Compruebo que el archivo existe.

                print(f"No se ha encontrado el archivo {ruta_archivo}")
                continue

            # Creo un conjunto para guardar los identificadores de las variantes que ya aparecen en la tabla de entrada.

            variantes_spliceai = set()

            # Me quedo con las variantes de la tabla generada previamente porque son las que ya han pasado los filtros anteriores.

            with open(tabla_entrada, "r", encoding="utf-8") as f_in:

                # La primera línea es la cabecera, así que la salto.

                next(f_in)

                for linea in f_in:

                    variante_id = linea.split("\t")[0]

                    variantes_spliceai.add(variante_id)

                # Abro el VCF anotado con VEP y recorro todas sus variantes.

                with gzip.open(ruta_archivo, "rt", encoding="utf-8") as f_vcf:

                    for variant in f_vcf:

                        if variant.startswith(
                        "#"
                        ):  # Las líneas que empiezan por # pertenecen a la cabecera y no las necesito.

                            continue

                        if not variant:  # Si la línea está vacía, la salto.

                            continue

                        cols = variant.split("\t")

                        # Formo el identificador de la variante para poder relacionarla con la tabla de entrada.

                        variante_id = f"{cols[0]}:{cols[1]}:{cols[3]}:{cols[4]}"

                        datos_info = utils.parse_info(cols[7])

                        csq_anot = datos_info.get("CSQ", ".")

                        # Obtengo los scores ADA y RF de dbscSNV.

                        ada_val, rf_val = obtener_scores_dbscsnv(csq_anot)

                        # Guardo los dos scores asociados al identificador de la variante.

                        memoria_dbscsnv[variante_id] = {
                            'ada': ada_val,
                            'rf': rf_val
                        }

                # Abro la tabla de entrada y creo la nueva tabla donde añadiré las dos clasificaciones de dbscSNV.

                with open(tabla_entrada, "r", encoding="utf-8") as f_in, open(tabla_salida, "w", encoding="utf-8") as f_out:

                    cabecera_original = f_in.readline().strip()  # Leo la cabecera y le quito el salto de línea.

                    # Añado las dos columnas nuevas a la cabecera.

                    f_out.write(
                        f"{cabecera_original}\tADA_dbscSNV(0.6)\tRF_dbscSNV(0.6)\n"
                    )

                    for linea in f_in:

                        linea_limpia = linea.strip()  # Quito el salto de línea.

                        columnas_tsv = linea_limpia.split("\t")  # Separo las columnas por tabuladores.

                        variante_id = columnas_tsv[0]  # La primera columna contiene la variante. Ejemplo: 14:55466853:A:G

                        # Busco la variante en el diccionario donde he guardado los scores.
                        # Si no la encuentro, utilizo -1 para indicar que no está anotada.

                        scores = memoria_dbscsnv.get(
                            variante_id,
                            {'ada': -1.0, 'rf': -1.0}
                        )

                        # Clasifico ADA como patogénica o benigna utilizando el umbral 0.6.

                        if scores["ada"] == -1.0:

                            class_ada = "No anotada"

                        else:

                            if scores["ada"] >= 0.6:

                                class_ada = "Patogénica"

                            else:

                                class_ada = "Benigna"

                        # Clasifico RF de la misma forma utilizando el umbral 0.6.

                        if scores["rf"] == -1.0:

                            class_rf = "No anotada"

                        else:

                            if scores["rf"] >= 0.6:

                                class_rf = "Patogénica"

                            else:

                                class_rf = "Benigna"

                        # Escribo la línea original y añado las dos clasificaciones de dbscSNV.

                        f_out.write(
                            f"{linea_limpia}\t{class_ada}\t{class_rf}\n"
                        )


if __name__ == "__main__":
    main()