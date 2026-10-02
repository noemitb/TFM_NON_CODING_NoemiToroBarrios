#!/usr/bin/env python3

# Noemí Toro Barrios | TFM: Revisión de herramientas y recursos para la anotación y priorización de variantes en regiones no codificantes.

# Script para extraer la anotación de SpliceAI, clasificar las variantes según los umbrales 0.2 y 0.8 y generar la tabla de resultados.

# Cargo todas las librerías, incluido el archivo donde tengo definidas las funciones pero sin el .py

import sys
import gzip
import os
import utils

# Nombre del archivo donde guardaré la tabla final

output = "megatabla_resultados_herramientas.tsv"


# Parseo la salida de SpliceAI y obtengo el score máximo de los cuatro valores.
# Si la variante no tiene anotación, le pongo -1 para poder detectarla después como No anotada.

def obtener_score_max(spliceai_str):

    if spliceai_str == ".":

        return -1

    try:  # Uso try por si hay algún error en la anotación. En principio no haría falta, pero lo dejo por si acaso.

        scores = []  # Aquí voy guardando los scores válidos.

        anotaciones = spliceai_str.split(",")

        for anotacion in anotaciones:

            partes = anotacion.split("|")  # Separo cada campo utilizando | como separador.

            if len(partes) < 6:  # Si faltan campos, considero que esta anotación no es válida.

                continue

            for valor in partes[
                2:6
            ]:  # Recorro los cuatro scores de SpliceAI. En Python el 6 no se incluye.

                if valor == "" or valor == ".":
                    continue  # Si el valor está vacío, lo salto.

                else:

                    scores.append(float(valor))

        if not scores:  # Si no he podido recuperar ningún score válido, devuelvo -1.

            return -1
        else:
            return max(scores)  # Si hay varios scores, me quedo con el valor máximo.

    except:

        return -1


# Compruebo que he lanzado el script indicando el archivo con la lista de VCFs.

def main():

    if len(sys.argv) < 2:

        print("No has indicado lista de VCFs para extraer la anotación")
        sys.exit(1)

    # Cargo el archivo que contiene las rutas de los VCFs que quiero procesar.

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

    # Creo el archivo donde voy a guardar todos los resultados.

    with open(output, "w", encoding="utf-8") as f_out:

        # Escribo la cabecera de la tabla.

        cabecera = "Variante\tConsecuencia\tClasificación Clinvar\tClasificación SpliceAI(0.2)\tClasificación SpliceAI(0.8)\n"
        f_out.write(cabecera)

        # Recorro uno a uno los archivos que tengo en la lista.

        for indice, ruta_archivo in enumerate(lista_rutas, start=1):

            nombre_archivo = os.path.basename(
                ruta_archivo
            )  # Me quedo solamente con el nombre del archivo para mostrarlo por pantalla.

            print(f"Leyendo y extrayendo datos de: {nombre_archivo}")

            if not os.path.exists(ruta_archivo):  # Compruebo que el archivo existe.

                print(f"No se ha encontrado el archivo {ruta_archivo}")
                continue

            with open(ruta_archivo, "r", encoding="utf-8") as f_vcf:  # Abro el VCF.

                for variant in f_vcf:

                    if variant.startswith(
                        "#"
                    ):  # Las líneas que empiezan por # pertenecen a la cabecera y no las necesito.

                        continue

                    if not variant:  # Si la línea está vacía, la salto.

                        continue

                    cols = variant.split("\t")  # Separo las columnas del VCF utilizando tabuladores.

                    info_raw = cols[7]  # Me quedo con la columna INFO.

                    datos_info = utils.parse_info(
                        info_raw
                    )  # Parseo la columna INFO utilizando la función que tengo definida en utils.

                    revstat = datos_info.get(
                        "CLNREVSTAT", "."
                    )  # Si no encuentro CLNREVSTAT, pongo "." para que las funciones puedan seguir trabajando.

                    clnsig = datos_info.get("CLNSIG")

                    clnsig_conf = datos_info.get("CLNSIGCONF")

                    # Clasifico la variante en uno de mis grupos de ClinVar.
                    # Si es Conflicting, también compruebo si contiene alguna clasificación patogénica.

                    categoria_clinvar, conflicting_pato = utils.clasificar_variante(
                        clnsig, clnsig_conf
                    )

                    # Obtengo el número de estrellas de ClinVar a partir del estado de revisión.

                    estrellas = utils.calcular_estrellas(revstat)

                    # Aplico el filtrado de ClinVar definido en utils.

                    if not utils.filtrar(
                        estrellas, conflicting_pato
                    ):  # Si la función devuelve False, la variante no pasa mis criterios y no la guardo.

                        continue

                    # Extraigo la anotación de SpliceAI. Si no existe, utilizo ".".

                    spliceai_anot = datos_info.get(
                        "SpliceAI", "."
                    )

                    # Obtengo el score máximo entre los cuatro scores de SpliceAI.

                    score_max = obtener_score_max(
                        spliceai_anot
                    )

                    # Clasifico cada variante utilizando los dos umbrales de SpliceAI que quiero estudiar: 0.2 y 0.8.

                    clasificacion_02 = ""

                    clasificacion_08 = ""

                    if score_max == -1:

                        clasificacion_02 = "No anotada"

                        clasificacion_08 = "No anotada"

                    else:

                        # Clasificación utilizando el umbral 0.2.

                        if score_max >= 0.2:

                            clasificacion_02 = "Patogénica"

                        else:

                            clasificacion_02 = "Benigna"

                        # Clasificación utilizando el umbral 0.8.

                        if score_max >= 0.8:

                            clasificacion_08 = "Patogénica"

                        else:

                            clasificacion_08 = "Benigna"

                    # Creo un identificador único de variante utilizando cromosoma, posición, referencia y alternativo.

                    variante_id = f"{cols[0]}:{cols[1]}:{cols[3]}:{cols[4]}"  # Ejemplo: 14:55466853:A:G

                    # Extraigo la consecuencia molecular de ClinVar.

                    mc_raw = datos_info.get(
                        "MC", "."
                    )

                    # Limpio las consecuencias para quedarme solamente con sus nombres y eliminar duplicados.

                    if mc_raw == "." or mc_raw == "NO HAY VALOR":

                        c_limpia = "Sin_Dato"

                    else:

                        consecuencias = mc_raw.split(",")

                        consecuencias_limpias = []

                        for c in consecuencias:

                            c = c.strip()

                            if "|" in c:

                                consecuencia = c.split("|")[-1]

                            else:

                                consecuencia = c

                            # Algunas entradas pueden contener varias consecuencias separadas por &.

                            for subcons in consecuencia.split("&"):

                                subcons = subcons.strip()

                                # Solo guardo la consecuencia si tiene contenido y todavía no la había añadido.

                                if subcons and subcons not in consecuencias_limpias:

                                    consecuencias_limpias.append(subcons)

                        if consecuencias_limpias:

                            # Si una variante tiene varias consecuencias, las guardo separadas por |.

                            c_limpia = "|".join(consecuencias_limpias)

                        else:

                            c_limpia = "Sin_Dato"

                    # Monto la línea final con toda la información que quiero guardar en la tabla.

                    linea_tsv = f"{variante_id}\t{c_limpia}\t{categoria_clinvar}\t{clasificacion_02}\t{clasificacion_08}\n"

                    f_out.write(linea_tsv)


if __name__ == "__main__":
    main()