#!/usr/bin/env python3

# Noemí Toro Barrios | TFM - Revisión de herramientas y recursos para la anotación y priorización de variantes en regiones no codificantes
#Clasificación de variantes según los niveles de GREEN-VARAN

#Cargo las librerías y las funciones definidas en utils.py

import sys
import gzip
import os
import utils

tabla_entrada = "megatabla_regulacion_fathmm_utrannotation.tsv"
tabla_salida = "megatabla_regulacion_completa.tsv"

#Establezco el umbral de nivel que voy a utilizar.

umbral_nivel = 3

#Obtengo la anotación de GREEN-VARAN del campo INFO.


def obtener_anotacion_greendb(datos_info):
    greendb_id = datos_info.get("greendb_id", ".")

    #La herramienta solo anota las variantes con las que solapa.

    if not greendb_id or greendb_id == "." or greendb_id == "NO HAY VALOR":
        return -1

    greendb_level_raw = datos_info.get("greendb_level", ".")
    
    try:
        return int(greendb_level_raw)
    except (TypeError, ValueError):
        return -1


def main():

    if len(sys.argv) < 2:
        print("Falta lista de vcfs para extraer la anotación")
        return

    archivo_lista = sys.argv[1]
    lista_rutas = []

    with open(archivo_lista, "r") as f_lista:  #Abro el archivo en modo lectura.
        for linea in f_lista:
            ruta = (
                linea.strip()
            )  #Quito saltos de línea o espacios.

            if ruta:
                lista_rutas.append(ruta)  #Añado la ruta a la lista si no está vacía.

    #Creo un diccionario para guardar los valores.

    memoria_greendb = {}

    for indice, ruta_archivo in enumerate(lista_rutas, start=1):

        #Me quedo con el nombre del archivo.
        nombre_archivo = os.path.basename(ruta_archivo)

        print(f"Leyendo y extrayendo datos de: {nombre_archivo}")

        if not os.path.exists(ruta_archivo):  #Compruebo que existe el archivo.
            print(f"No se ha encontrado el archivo {ruta_archivo}")
            continue

        with gzip.open(ruta_archivo, "rt", encoding="utf-8") as f_vcf:  #Abro el VCF.
            for variant in f_vcf:
                #Si comienza con almohadilla no me interesa porque es la cabecera.
                if variant.startswith("#"):
                    continue

                if not variant:  #Si está vacía, la ignoro.
                    continue

                cols = variant.split("\t")
                cromosoma = cols[0].replace("chr", "") if cols[0].startswith("chr") else cols[0]
                variante_id = f"{cromosoma}:{cols[1]}:{cols[3]}:{cols[4]}"  #Formo el identificador de la variante.
                datos_info = utils.parse_info(cols[7])

                nivel_val = obtener_anotacion_greendb(datos_info)
                if variante_id not in memoria_greendb or memoria_greendb[variante_id] == -1:
                    memoria_greendb[variante_id] = nivel_val
                elif nivel_val != -1:
                    memoria_greendb[variante_id] = max(memoria_greendb[variante_id], nivel_val)
                
    #Abro la tabla generada con las otras herramientas.
    
    with open(tabla_entrada, "r", encoding="utf-8") as f_in, open(tabla_salida, "w", encoding="utf-8") as f_out:

        cabecera_original = f_in.readline().strip()  #Leo la primera línea y le quito el salto de línea.

        #Añado la columna de GREEN-VARAN a la cabecera.
        nombre_columna = f"GREEN-DB({umbral_nivel})"
        f_out.write(f"{cabecera_original}\t{nombre_columna}\n")

        for linea in f_in:

            linea_limpia = linea.strip()  #Quito el salto de línea.
            columnas_tsv = linea_limpia.split("\t")  #Separo por tabuladores.
            variante_id = columnas_tsv[0]  #La primera columna es la variante.

            #Busco esta variante en el diccionario creado con los niveles de GREEN-VARAN.
            nivel = memoria_greendb.get(variante_id, -1)  #Si no está, indico que no está anotada.

            #Clasifico las variantes utilizando el umbral de GREEN-VARAN.
            if nivel == -1:
                clasificacion = "No anotada"
            elif nivel >= umbral_nivel:
                clasificacion = "Patogénica"
            else:
                clasificacion = "Benigna"

            #Escribo la línea original junto con la columna de GREEN-VARAN.
            f_out.write(f"{linea_limpia}\t{clasificacion}\n")


if __name__ == "__main__":
    main()