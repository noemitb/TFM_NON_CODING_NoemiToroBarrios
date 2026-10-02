#!/usr/bin/env python3

# Noemí Toro Barrios | TFM - Revisión de herramientas y recursos para la anotación y priorización de variantes en regiones no codificantes
#Generación de la tabla de FATHMM-MKL y UTRAnnotator

import sys
import gzip
import os
import re
import utils


output = "megatabla_regulacion_fathmm_utrannotation.tsv"


def abrir_vcf(ruta):

    if ruta.endswith(".gz"):
        return gzip.open(ruta, "rt", encoding="utf-8")

    return open(ruta, "r", encoding="utf-8")


def leer_lista(ruta_lista):

    lista_rutas = []

    with open(ruta_lista, "r") as f_lista:

        for linea in f_lista:

            ruta = linea.strip()

            if ruta:
                lista_rutas.append(ruta)

    return lista_rutas


def normalizar_chr(chrom):

    chrom = str(chrom)

    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]

    return chrom


def crear_clave(chrom, pos, ref, alt):

    return (
        normalizar_chr(chrom),
        str(pos),
        ref.upper(),
        alt.upper()
    )


def obtener_consecuencias(mc_raw):

    #Si no hay consecuencia molecular, indico que no hay dato.
    if not mc_raw or mc_raw == "." or mc_raw == "NO HAY VALOR":

        return "Sin_Dato"

    consecuencias = []

    #Separo las distintas consecuencias de ClinVar.
    for c in mc_raw.split(","):

        c = c.strip()

        if not c:
            continue

        #El campo MC tiene el formato:
        #SO:0001583|missense_variant
        if "|" in c:

            c_limpia = c.split("|")[-1]

        else:

            c_limpia = c

        #Si aparecen varias consecuencias separadas por &, las separo.
        for subcons in c_limpia.split("&"):

            subcons = subcons.strip()

            if subcons and subcons not in consecuencias:

                consecuencias.append(subcons)

    if not consecuencias:

        return "Sin_Dato"

    #Mantengo todas las consecuencias en la misma columna.
    return "|".join(consecuencias)


def obtener_score_fathmm(valor):

    if not valor or valor == ".":
        return None

    scores = []

    #Si hay más de un valor, me quedo con el máximo.
    for score in re.split(r"[,|&]", valor):

        try:

            scores.append(float(score))

        except ValueError:

            continue

    if scores:

        return max(scores)

    return None


def cargar_fathmm(lista_green):

    memoria_fathmm = {}

    lista_rutas = leer_lista(lista_green)

    total = 0
    anotadas = 0

    for ruta_archivo in lista_rutas:

        nombre_archivo = os.path.basename(ruta_archivo)

        if not os.path.exists(ruta_archivo):

            continue

        with abrir_vcf(ruta_archivo) as f_vcf:

            for variant in f_vcf:

                if variant.startswith("#"):
                    continue

                if not variant:
                    continue

                cols = variant.rstrip("\n").split("\t")

                if len(cols) < 8:
                    continue

                chrom = cols[0]
                pos = cols[1]
                ref = cols[3]
                alt = cols[4]

                #Descarto las variantes multialélicas.
                if "," in alt:
                    continue

                total += 1

                datos_info = utils.parse_info(cols[7])

                fathmm_raw = datos_info.get(
                    "FATHMM_MKLNC",
                    "."
                )

                score = obtener_score_fathmm(
                    fathmm_raw
                )

                if score is None:
                    continue

                clave = crear_clave(
                    chrom,
                    pos,
                    ref,
                    alt
                )

                memoria_fathmm[clave] = score

                anotadas += 1

    return memoria_fathmm


def obtener_loeuf(datos_info, info_raw):

    loeuf_raw = datos_info.get(
        "LOEUF_score"
    )

    if loeuf_raw and loeuf_raw != ".":

        try:

            return float(loeuf_raw)

        except ValueError:

            pass

    #Busco LOEUF directamente en INFO por seguridad.
    match = re.search(
        r"(?:^|;)LOEUF_score=([^;]+)",
        info_raw
    )

    if match:

        try:

            return float(match.group(1))

        except ValueError:

            pass

    return None


def obtener_campos_csq(linea):

    match = re.search(
        r"Format:\s*([^\">]+)",
        linea
    )

    if not match:

        return None

    return (
        match
        .group(1)
        .strip()
        .split("|")
    )


def obtener_datos_utr_completos(csq_str, campos_csq):

    if not csq_str or csq_str == ".":
        return None, None

    if not campos_csq:
        return None, None

    try:

        indice_annot = campos_csq.index(
            "5UTR_annotation"
        )

        indice_cons = campos_csq.index(
            "5UTR_consequence"
        )

    except ValueError:

        return None, None

    utr_annot_raw = None
    utr_cons_raw = None

    #Si hay varios transcritos con anotación válida, me quedo con la última anotación UTR válida.

    for transcrito in csq_str.split(","):

        partes = transcrito.split("|")

        if len(partes) <= max(
            indice_annot,
            indice_cons
        ):

            continue

        utrannotation_csq = (
            partes[indice_annot].strip()
        )

        utrannotation_tag = (
            partes[indice_cons].strip()
        )

        if (
            utrannotation_tag
            and utrannotation_tag != "."
        ):

            utr_annot_raw = utrannotation_csq
            utr_cons_raw = utrannotation_tag

    return utr_annot_raw, utr_cons_raw


def clasificar_utr_zhang(
    utr_annot,
    utr_cons,
    loeuf_val
):

    if (
        loeuf_val is None
        or not utr_annot
        or utr_annot == "."
    ):

        return "No anotado"

    annot_lower = utr_annot.lower()

    cons_lower = (
        str(utr_cons).lower()
        if utr_cons
        else ""
    )

    is_overlapping = (
        "oorf" in cons_lower
        or "oorf" in annot_lower
    )

    texto_total = (
        annot_lower
        + " "
        + cons_lower
    )

    kozak = "Unknown"

    if "strong" in texto_total:

        kozak = "Strong"

    elif "moderate" in texto_total:

        kozak = "Moderate"

    elif "weak" in texto_total:

        kozak = "Weak"

    if (
        is_overlapping
        and kozak in ["Strong", "Moderate"]
        and loeuf_val < 0.35
    ):

        return "patogénico"

    elif (
        is_overlapping
        and kozak == "Weak"
        and loeuf_val < 0.35
    ):

        return "vus"

    elif (
        not is_overlapping
        and kozak == "Strong"
        and loeuf_val < 0.35
    ):

        return "probablemente patogénico"

    elif (
        is_overlapping
        and kozak == "Strong"
        and loeuf_val > 1.0
    ):

        return "probablemente benigno"

    else:

        return "benigno"


def main():

    if len(sys.argv) < 3:

        sys.exit(1)

    lista_utr = sys.argv[1]
    lista_green = sys.argv[2]

    #Cargo en memoria los scores de FATHMM-MKL de GREEN-VARAN.

    memoria_fathmm = cargar_fathmm(
        lista_green
    )

    rutas_utr = leer_lista(
        lista_utr
    )

    total_variantes = 0

    total_snv = 0
    total_indel = 0

    fathmm_anotadas = 0
    fathmm_no_anotadas = 0

    fathmm_snv = 0
    fathmm_indel = 0

    loeuf_presente = 0
    loeuf_ausente = 0

    utr_anotadas = 0
    utr_no_anotadas = 0

    multiples_consecuencias = 0

    with open(
        output,
        "w",
        encoding="utf-8"
    ) as f_out:

        #Mantengo la estructura original de la tabla.
        f_out.write(
            "Variante\t"
            "Consecuencia\t"
            "Clasificación Clinvar\t"
            "FATHH\t"
            "UTRannotator\n"
        )

        for ruta_archivo in rutas_utr:

            nombre_archivo = os.path.basename(
                ruta_archivo
            )

            if not os.path.exists(
                ruta_archivo
            ):

                continue

            campos_csq = None

            with abrir_vcf(
                ruta_archivo
            ) as f_vcf:

                for variant in f_vcf:

                    if variant.startswith(
                        "##INFO=<ID=CSQ"
                    ):

                        campos_csq = (
                            obtener_campos_csq(
                                variant
                            )
                        )

                        continue

                    if variant.startswith("#"):
                        continue

                    if not variant:
                        continue

                    cols = (
                        variant
                        .rstrip("\n")
                        .split("\t")
                    )

                    if len(cols) < 8:
                        continue

                    chrom = cols[0]
                    pos = cols[1]
                    ref = cols[3]
                    alt = cols[4]
                    info_raw = cols[7]

                    if "," in alt:

                        continue

                    datos_info = (
                        utils.parse_info(
                            info_raw
                        )
                    )

                    clnsig = datos_info.get(
                        "CLNSIG",
                        ""
                    )

                    clnsig_conf = (
                        datos_info.get(
                            "CLNSIGCONF",
                            ""
                        )
                    )

                    revstat = datos_info.get(
                        "CLNREVSTAT",
                        ""
                    )

                    clasificacion, patogenico = (
                        utils.clasificar_variante(
                            clnsig,
                            clnsig_conf
                        )
                    )

                    estrellas = (
                        utils.calcular_estrellas(
                            revstat
                        )
                    )

                    if not utils.filtrar(
                        estrellas,
                        patogenico
                    ):

                        continue

                    total_variantes += 1

                    mc_raw = datos_info.get(
                        "MC",
                        "."
                    )

                    consecuencia = (
                        obtener_consecuencias(
                            mc_raw
                        )
                    )

                    if "|" in consecuencia:

                        multiples_consecuencias += 1

                    if (
                        len(ref) == 1
                        and len(alt) == 1
                    ):

                        total_snv += 1
                        tipo = "SNV"

                    else:

                        total_indel += 1
                        tipo = "INDEL"

                    clave = crear_clave(
                        chrom,
                        pos,
                        ref,
                        alt
                    )

                    fathmm_score = (
                        memoria_fathmm.get(
                            clave
                        )
                    )

                    if fathmm_score is None:

                        class_fathmm = (
                            "No anotado"
                        )

                        fathmm_no_anotadas += 1

                    else:

                        fathmm_anotadas += 1

                        if tipo == "SNV":

                            fathmm_snv += 1

                        else:

                            fathmm_indel += 1

                        if fathmm_score >= 0.5:

                            class_fathmm = (
                                "patogénico"
                            )

                        else:

                            class_fathmm = (
                                "benigno"
                            )

                    loeuf_val = obtener_loeuf(
                        datos_info,
                        info_raw
                    )

                    if loeuf_val is None:

                        loeuf_ausente += 1

                    else:

                        loeuf_presente += 1

                    csq_str = datos_info.get(
                        "CSQ",
                        "."
                    )

                    (
                        utr_annot,
                        utr_cons
                    ) = obtener_datos_utr_completos(
                        csq_str,
                        campos_csq
                    )

                    class_utr = (
                        clasificar_utr_zhang(
                            utr_annot,
                            utr_cons,
                            loeuf_val
                        )
                    )

                    if (
                        class_utr
                        == "No anotado"
                    ):

                        utr_no_anotadas += 1

                    else:

                        utr_anotadas += 1

                    variante_id = (
                        f"{chrom}:"
                        f"{pos}:"
                        f"{ref}:"
                        f"{alt}"
                    )

                    f_out.write(
                        f"{variante_id}\t"
                        f"{consecuencia}\t"
                        f"{clasificacion}\t"
                        f"{class_fathmm}\t"
                        f"{class_utr}\n"
                    )


if __name__ == "__main__":

    main()