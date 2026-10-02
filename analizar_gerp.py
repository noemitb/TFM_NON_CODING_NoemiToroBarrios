#!/usr/bin/env python3

# Noemí Toro Barrios | TFM - Revisión de herramientas y recursos para la anotación y priorización de variantes en regiones no codificantes
#Clasificación de variantes según las puntuaciones de GERP++

#Cargo las librerías y las funciones definidas en utils.py

import sys
import gzip
import os
import utils

tabla_entrada = "megatabla_regulación_completa_cadd.tsv"
tabla_salida = "megatabla_regulación_completa_cadd_gerp.tsv"

#Parseo la salida de VEP para obtener el resultado de GERP++.

def obtener_scores_gerp(csq_str): #La anotación de la herramienta aparece en el campo CSQ

    if not csq_str or csq_str == ".": #Compruebo que esté anotado ese campo.
        return -10000.0 #Indico que la variante no está anotada
    
    
    transcritos = csq_str.split(',') #Como puede haber diferentes anotaciones por transcrito, separo cada uno de los transcritos.
    scores_gerp = []

    for t in transcritos:
        partes = t.split('|')

        if len(partes) <= 32:
            continue

        gerp_raw = partes[32].strip()
        

        if gerp_raw and gerp_raw != "." and gerp_raw != "": #Convierto el valor en un número.

            try:
                scores_gerp.append(float(gerp_raw))

            except ValueError:
                continue
    
            
    #Me quedo con el valor máximo del predictor si es diferente de -

    if scores_gerp:
        return max(scores_gerp) 
    else:
        return -10000.0


def main():

    if len(sys.argv) < 2:

        print("No has indicado lista de VCFs para extraer la anotación")
        sys.exit(1)

    #Si hay input, comienzo con la generación de tablas.
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

    memoria_gerp = {}

    for indice, ruta_archivo in enumerate(lista_rutas, start=1):

            nombre_archivo = os.path.basename(ruta_archivo)  #Me quedo con el nombre del archivo.

            print(f"Leyendo y extrayendo datos de: {nombre_archivo}")

            if not os.path.exists(ruta_archivo):  #Compruebo que existe el archivo.

                print(f"No se ha encontrado el archivo {ruta_archivo}")
                continue

            variantes_spliceai = set()

#Me quedo con los identificadores de las variantes de la tabla de entrada.
            with open(tabla_entrada,"r", encoding="utf-8") as f_in:
                #La primera línea es la cabecera así que la salto.
                next(f_in)

                for linea in f_in:
                    variante_id = linea.split("\t")[0]
                    variantes_spliceai.add(variante_id)

            with gzip.open(ruta_archivo, "rt", encoding="utf-8") as f_vcf:  #Abro el VCF.
                
                    for variant in f_vcf:

                        if variant.startswith(
                        "#"
                        ):  #Si comienza con almohadilla no me interesa porque es la cabecera.
                            continue

                        if not variant:  #Si está vacía, la ignoro.
                            continue

                        cols = variant.rstrip("\n").split("\t")

                        if len(cols) < 8:
                            continue

                        variante_id = f"{cols[0]}:{cols[1]}:{cols[3]}:{cols[4]}" #Formo el identificador de la variante.

                        datos_info = utils.parse_info(cols[7])

                        csq_anot = datos_info.get("CSQ", ".")

                        gerp_val = obtener_scores_gerp(csq_anot)

                        memoria_gerp[variante_id] = gerp_val

    #Abro la tabla de entrada y la nueva tabla de salida.
    with open(tabla_entrada, "r", encoding="utf-8") as f_in, open(tabla_salida, "w", encoding="utf-8") as f_out:
        
        cabecera_original = f_in.readline().strip() #Leo la primera línea y le quito el salto de línea.
        
        #Añado la columna de GERP++ a la cabecera.

        f_out.write(f"{cabecera_original}\tGERP(2)\n")
        
        for linea in f_in:
                        
            linea_limpia = linea.strip() #Quito el salto de línea.

            columnas_tsv = linea_limpia.split("\t") #Separo por tabuladores.

            variante_id = columnas_tsv[0] #La primera columna es la variante.
            
            #Busco esta variante en el diccionario creado con los scores de la herramienta.
            scores = memoria_gerp.get(variante_id,-10000.0) #Si no está, indico que no está anotada.
            
            #Clasifico las variantes utilizando el umbral de GERP++.

            if scores == -10000.0: 
                class_gerp = "No anotada"

            else: 

                if scores >= 2:
                    class_gerp = "Patogénica" 

                else:
                    class_gerp = "Benigna"
                            
                            
            #Escribo la línea original junto con la columna de GERP++.
            f_out.write(f"{linea_limpia}\t{class_gerp}\n")


if __name__ == "__main__":
    main()