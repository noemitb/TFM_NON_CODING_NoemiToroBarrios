#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd

# Noemí Toro Barrios | TFM: Revisión de herramientas y recursos para la anotación y priorización de variantes en regiones no codificantes

#Script para realizar el análisis comparativo. 


#Defino los umbrales de CADD y SpliceAI

#Guardo las columnas de cada umbral
CADD_CONFIG = {
    10: "CADD(10)",
    15: "CADD(15)",
    20: "CADD(20)",
}

SPLICEAI_CONFIG = {
    0.2: "Clasificación SpliceAI(0.2)",
    0.8: "Clasificación SpliceAI(0.8)",
}

DBSCSNV_MODE = "BOTH"
STRICT_NONCODING_FILTER = True

#Defino las consecuencias codificantes que voy a excluir
CODING_CONSEQUENCES = {
    "missense_variant",
    "synonymous_variant",
    "frameshift_variant",
    "nonsense",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "initiator_codon_variant",
    "inframe_insertion",
    "inframe_deletion",
    "inframe_indel",
    "protein_altering_variant",
    "coding_sequence_variant",
    "incomplete_terminal_codon_variant",
}

#Defino las consecuencias que voy a comparar
TARGET_CONSEQUENCES = [
    "3_prime_UTR_variant",
    "5_prime_UTR_variant",
    "intron_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "genic_upstream_transcript_variant",
    "genic_downstream_transcript_variant",
]

Z_975 = 1.959963984540054


#Creo las funciones auxiliares

def normalizar_texto(valor):
    if pd.isna(valor):
        return ""
    s = str(valor).strip().lower()
    return "".join(
        c for c in unicodedata.normalize("NFD", s)
        if unicodedata.category(c) != "Mn"
    )


def binarizar_prediccion(valor):
    #Paso las predicciones a P o B
    #Si no hay clasificación válida dejo NA
    s = normalizar_texto(valor)

    if not s or s in {"nan", "na", "n/a", ".", "-"}:
        return pd.NA

    if "patogen" in s or "pathogen" in s:
        return "P"

    if "benign" in s:
        return "B"

    return pd.NA


def categorizar_clinvar(valor):
    s = normalizar_texto(valor)

    if "conflict" in s:
        return "Conflicting"

    if "vus" in s or "inciert" in s or "uncertain" in s:
        return "VUS"

    if "patogen" in s or "pathogen" in s:
        return "P/LP"

    if "benign" in s:
        return "B/LB"

    return "Otras"


def pct(numerador, denominador):
    if denominador == 0:
        return np.nan
    return round(100.0 * numerador / denominador, 2)


def separar_consecuencias(valor):
    if pd.isna(valor):
        return []
    return [x.strip() for x in str(valor).split("|") if x.strip()]


def herramientas_dbscsnv():
    mode = DBSCSNV_MODE.upper()

    if mode == "BOTH":
        return ["ADA_dbscSNV(0.6)", "RF_dbscSNV(0.6)"]
    if mode == "ADA":
        return ["ADA_dbscSNV(0.6)"]
    if mode == "RF":
        return ["RF_dbscSNV(0.6)"]

    raise ValueError("DBSCSNV_MODE debe ser 'BOTH', 'ADA' o 'RF'.")


#Calculo los intervalos de confianza

def wilson_interval(k, n, z=Z_975):
    #Calculo el IC95% de Wilson
    if n is None or n <= 0:
        return np.nan, np.nan

    k = int(k)
    n = int(n)

    if k < 0 or k > n:
        return np.nan, np.nan

    p = k / n
    z2 = z * z
    den = 1.0 + z2 / n

    centro = (p + z2 / (2.0 * n)) / den
    semi = (
        z
        * np.sqrt(
            p * (1.0 - p) / n
            + z2 / (4.0 * n * n)
        )
        / den
    )

    inf = max(0.0, centro - semi) * 100.0
    sup = min(1.0, centro + semi) * 100.0

    return round(inf, 2), round(sup, 2)


def bootstrap_ganancia_pareada(
    preds_p,
    n_bootstrap,
    rng,
):
    #Hago el bootstrap con las variantes P/LP
    #Recalculo la sensibilidad en cada réplica
    #Calculo la ganancia frente a la mejor herramienta
    #Agrupo patrones P/B para acelerar el cálculo
    if preds_p is None or len(preds_p) == 0:
        return np.nan, np.nan, np.nan

    n = len(preds_p)
    m = preds_p.shape[1]

    #Convierto P y B en 1 y 0
    arr = (preds_p.to_numpy() == "P").astype(np.int8)

    #Busco los patrones repetidos
    patterns, counts = np.unique(arr, axis=0, return_counts=True)
    probs = counts / counts.sum()

    #Compruebo si alguna herramienta predice P
    or_pattern = patterns.any(axis=1).astype(np.int8)

    #Hago el bootstrap por bloques
    chunk = min(2000, n_bootstrap)
    ganancias = np.empty(n_bootstrap, dtype=float)

    start = 0
    while start < n_bootstrap:
        b = min(chunk, n_bootstrap - start)

        #Genero los conteos bootstrap
        boot_counts = rng.multinomial(
            n=n,
            pvals=probs,
            size=b,
        )

        #Calculo la sensibilidad de cada herramienta
        det_ind = boot_counts @ patterns  #Cada fila es una réplica y cada columna una herramienta
        sens_ind = det_ind / n

        #Me quedo con la mejor sensibilidad
        best_ind = sens_ind.max(axis=1)

        #Calculo la sensibilidad de la estrategia de unión
        det_or = boot_counts @ or_pattern
        sens_or = det_or / n

        ganancias[start:start+b] = 100.0 * (sens_or - best_ind)
        start += b

    inf, mediana, sup = np.percentile(
        ganancias,
        [2.5, 50.0, 97.5],
    )

    return (
        round(float(inf), 2),
        round(float(sup), 2),
        round(float(mediana), 2),
    )


#Preparo los datos para el análisis

def preparar_datos(df):
    columnas_necesarias = {
        "Variante",
        "Consecuencia",
        "Clasificación Clinvar",
        "Tipo_variante",
        "FATHH",
        "UTRannotator",
        "GREEN-DB(3)",
        "CADD(10)",
        "CADD(15)",
        "CADD(20)",
        "GERP(2)",
        "Clasificación SpliceAI(0.2)",
        "Clasificación SpliceAI(0.8)",
        "ADA_dbscSNV(0.6)",
        "RF_dbscSNV(0.6)",
        "PDIVAS(0.082)",
    }

    faltan = sorted(columnas_necesarias - set(df.columns))
    if faltan:
        raise ValueError(
            "Faltan columnas necesarias en tabla_fusionada.tsv:\n  - "
            + "\n  - ".join(faltan)
        )

    df = df.copy()

    #Separo las consecuencias anotadas
    df["_lista_consecuencias"] = (
        df["Consecuencia"].apply(separar_consecuencias)
    )

    #Excluyo variantes con alguna consecuencia codificante
    if STRICT_NONCODING_FILTER:
        es_coding = df["_lista_consecuencias"].apply(
            lambda xs: any(x in CODING_CONSEQUENCES for x in xs)
        )

        variantes_coding = set(
            df.loc[es_coding, "Variante"]
            .dropna()
            .astype(str)
        )

        if variantes_coding:
            df = df.loc[
                ~df["Variante"].astype(str).isin(variantes_coding)
            ].copy()

    #Dejo una fila por variante y consecuencia
    df = df.explode("_lista_consecuencias")
    df["Consecuencia"] = df["_lista_consecuencias"]
    df = df.drop(columns="_lista_consecuencias")

    df = df[
        df["Consecuencia"].isin(TARGET_CONSEQUENCES)
    ].copy()

    df["Tipo_variante"] = (
        df["Tipo_variante"]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    #Elimino duplicados
    df = df.drop_duplicates(
        subset=["Variante", "Consecuencia", "Tipo_variante"],
        keep="first",
    )

    df["_ClinVar_categoria"] = (
        df["Clasificación Clinvar"].apply(categorizar_clinvar)
    )

    pred_cols = [
        "FATHH",
        "UTRannotator",
        "GREEN-DB(3)",
        "CADD(10)",
        "CADD(15)",
        "CADD(20)",
        "GERP(2)",
        "Clasificación SpliceAI(0.2)",
        "Clasificación SpliceAI(0.8)",
        "ADA_dbscSNV(0.6)",
        "RF_dbscSNV(0.6)",
        "PDIVAS(0.082)",
    ]

    for col in pred_cols:
        df[f"_BIN__{col}"] = (
            df[col].apply(binarizar_prediccion)
        )

    return df


#Defino las herramientas de cada comparación

def herramientas_regulacion(consecuencia, tipo, cadd_col):
    if tipo == "SNV":
        if consecuencia == "5_prime_UTR_variant":
            return [
                "FATHH",
                "UTRannotator",
                "GREEN-DB(3)",
                cadd_col,
                "GERP(2)",
            ]

        return [
            "FATHH",
            "GREEN-DB(3)",
            cadd_col,
            "GERP(2)",
        ]

    return [
        "GREEN-DB(3)",
        cadd_col,
        "GERP(2)",
    ]


def herramientas_splicing(
    consecuencia,
    tipo,
    cadd_col,
    spliceai_col,
):
    herramientas = [spliceai_col]

    if consecuencia == "intron_variant":
        herramientas += ["PDIVAS(0.082)"]

    elif (
        consecuencia in {
            "splice_acceptor_variant",
            "splice_donor_variant",
        }
        and tipo == "SNV"
    ):
        herramientas += herramientas_dbscsnv()

    herramientas += [cadd_col, "GERP(2)"]
    return herramientas


#Analizo cada consecuencia y tipo de variante

def analizar_estrato(
    df,
    consecuencia,
    tipo,
    herramientas,
    umbral_cadd,
    umbral_spliceai,
    n_bootstrap,
    rng,
):
    sub = df[
        (df["Consecuencia"] == consecuencia)
        & (df["Tipo_variante"] == tipo)
    ].copy()

    n_total = len(sub)
    bin_cols = [f"_BIN__{h}" for h in herramientas]

    if n_total:
        mask_validas = sub[bin_cols].notna().all(axis=1)
        evaluadas = sub.loc[mask_validas].copy()
    else:
        evaluadas = sub.copy()

    n_eval = len(evaluadas)

    #Inicializo los contadores
    n_consenso = 0
    n_consenso_p = 0
    n_consenso_b = 0
    n_consenso_correcto_clinvar = 0

    dist_p = {
        "P/LP": 0,
        "VUS": 0,
        "Conflicting": 0,
        "B/LB": 0,
        "Otras": 0,
    }

    n_p_lp = 0
    n_p_lp_detectadas_or = 0
    n_b_lb = 0
    n_b_lb_fp_or = 0

    sensibilidad_or = np.nan
    fp_or = np.nan
    especificidad_or = np.nan
    mejor_sens_ind = np.nan
    ganancia_or = np.nan
    sensibilidades_por_herramienta = {}

    #Inicializo los IC
    conc_h_inf = conc_h_sup = np.nan
    conc_cv_inf = conc_cv_sup = np.nan
    sens_inf = sens_sup = np.nan
    fp_inf = fp_sup = np.nan
    vus_inf = vus_sup = np.nan
    conf_inf = conf_sup = np.nan
    gain_inf = gain_sup = gain_mediana = np.nan

    if n_eval:
        preds = evaluadas[bin_cols]

        #Calculo el consenso entre herramientas
        todas_p = (preds == "P").all(axis=1)
        todas_b = (preds == "B").all(axis=1)
        consenso = todas_p | todas_b

        n_consenso_p = int(todas_p.sum())
        n_consenso_b = int(todas_b.sum())
        n_consenso = int(consenso.sum())

        #Comparo el consenso con ClinVar
        clin_cat = evaluadas["_ClinVar_categoria"]

        consenso_p_correcto = (
            todas_p & (clin_cat == "P/LP")
        )
        consenso_b_correcto = (
            todas_b & (clin_cat == "B/LB")
        )

        n_consenso_correcto_clinvar = int(
            (
                consenso_p_correcto
                | consenso_b_correcto
            ).sum()
        )

        #Reviso ClinVar en los consensos P
        clinvar_consenso_p = evaluadas.loc[
            todas_p,
            "_ClinVar_categoria",
        ]

        conteos_p = (
            clinvar_consenso_p
            .value_counts(dropna=False)
            .to_dict()
        )

        for categoria in dist_p:
            dist_p[categoria] = int(
                conteos_p.get(categoria, 0)
            )

        #Considero P si alguna herramienta predice P
        positiva_or = (preds == "P").any(axis=1)

        mask_p_lp = clin_cat == "P/LP"
        mask_b_lb = clin_cat == "B/LB"

        n_p_lp = int(mask_p_lp.sum())
        n_b_lb = int(mask_b_lb.sum())

        n_p_lp_detectadas_or = int(
            (positiva_or & mask_p_lp).sum()
        )
        n_b_lb_fp_or = int(
            (positiva_or & mask_b_lb).sum()
        )

        sensibilidad_or = pct(
            n_p_lp_detectadas_or,
            n_p_lp,
        )
        fp_or = pct(
            n_b_lb_fp_or,
            n_b_lb,
        )

        if not pd.isna(fp_or):
            especificidad_or = round(
                100.0 - fp_or,
                2,
            )

        #Calculo la sensibilidad de cada herramienta
        #Me quedo con la sensibilidad más alta
        sensibilidades_individuales = []

        if n_p_lp:
            for herramienta, bc in zip(
                herramientas,
                bin_cols,
            ):
                detectadas = int(
                    (
                        (preds[bc] == "P")
                        & mask_p_lp
                    ).sum()
                )
                sens = pct(detectadas, n_p_lp)
                if not pd.isna(sens):
                    sensibilidades_individuales.append(
                        sens
                    )
                    sensibilidades_por_herramienta[herramienta] = sens

        if sensibilidades_individuales:
            mejor_sens_ind = max(
                sensibilidades_individuales
            )

            if not pd.isna(sensibilidad_or):
                ganancia_or = round(
                    sensibilidad_or
                    - mejor_sens_ind,
                    2,
                )

        #Calculo los IC95% de Wilson
        conc_h_inf, conc_h_sup = wilson_interval(
            n_consenso,
            n_eval,
        )

        conc_cv_inf, conc_cv_sup = wilson_interval(
            n_consenso_correcto_clinvar,
            n_eval,
        )

        sens_inf, sens_sup = wilson_interval(
            n_p_lp_detectadas_or,
            n_p_lp,
        )

        fp_inf, fp_sup = wilson_interval(
            n_b_lb_fp_or,
            n_b_lb,
        )

        vus_inf, vus_sup = wilson_interval(
            dist_p["VUS"],
            n_consenso_p,
        )

        conf_inf, conf_sup = wilson_interval(
            dist_p["Conflicting"],
            n_consenso_p,
        )

        #Calculo el IC95% de la ganancia
        if n_p_lp > 0:
            preds_p = preds.loc[
                mask_p_lp,
                bin_cols,
            ].copy()

            (
                gain_inf,
                gain_sup,
                gain_mediana,
            ) = bootstrap_ganancia_pareada(
                preds_p=preds_p,
                n_bootstrap=n_bootstrap,
                rng=rng,
            )

    #Preparo las sensibilidades individuales
    #Uso "-" si no se evalúa y "NA" si no se puede calcular
    def sens_individual_salida(nombre_exacto=None, prefijo=None):
        if nombre_exacto is not None:
            candidatos = [
                h for h in herramientas
                if h == nombre_exacto
            ]
        else:
            candidatos = [
                h for h in herramientas
                if h.startswith(prefijo)
            ]

        if not candidatos:
            return "-"

        sens = sensibilidades_por_herramienta.get(
            candidatos[0],
            np.nan,
        )
        if pd.isna(sens):
            return "NA"

        return round(float(sens), 2)

    return {
        "Consecuencia": consecuencia,
        "Tipo": tipo,
        "Umbral_CADD": umbral_cadd,
        "Umbral_SpliceAI": (
            umbral_spliceai
            if umbral_spliceai is not None
            else "NA"
        ),
        "Herramientas_comparadas":
            " + ".join(herramientas),

        "N_total": n_total,
        "N_evaluadas": n_eval,
        "Cobertura_conjunta_pct":
            pct(n_eval, n_total),

        "N_consenso_herramientas":
            n_consenso,
        "Concordancia_herramientas_pct":
            pct(n_consenso, n_eval),
        "Concordancia_herramientas_IC95_inf_pct":
            conc_h_inf,
        "Concordancia_herramientas_IC95_sup_pct":
            conc_h_sup,

        "N_consenso_correcto_ClinVar":
            n_consenso_correcto_clinvar,
        "Concordancia_consenso_ClinVar_sobre_Nevaluadas_pct":
            pct(
                n_consenso_correcto_clinvar,
                n_eval,
            ),
        "Concordancia_ClinVar_IC95_inf_pct":
            conc_cv_inf,
        "Concordancia_ClinVar_IC95_sup_pct":
            conc_cv_sup,

        "N_consenso_patogenico":
            n_consenso_p,
        "Consenso_patogenico_sobre_Nevaluadas_pct":
            pct(n_consenso_p, n_eval),

        "Consenso_P_ClinVar_VUS_N":
            dist_p["VUS"],
        "Consenso_P_ClinVar_VUS_pct":
            pct(
                dist_p["VUS"],
                n_consenso_p,
            ),
        "Consenso_P_ClinVar_VUS_IC95_inf_pct":
            vus_inf,
        "Consenso_P_ClinVar_VUS_IC95_sup_pct":
            vus_sup,

        "Consenso_P_ClinVar_Conflicting_N":
            dist_p["Conflicting"],
        "Consenso_P_ClinVar_Conflicting_pct":
            pct(
                dist_p["Conflicting"],
                n_consenso_p,
            ),
        "Consenso_P_ClinVar_Conflicting_IC95_inf_pct":
            conf_inf,
        "Consenso_P_ClinVar_Conflicting_IC95_sup_pct":
            conf_sup,

        "N_P_LP_ClinVar":
            n_p_lp,
        "N_P_LP_detectadas_por_ge1":
            n_p_lp_detectadas_or,
        "Sensibilidad_OR_ge1_pct":
            sensibilidad_or,
        "Sensibilidad_OR_IC95_inf_pct":
            sens_inf,
        "Sensibilidad_OR_IC95_sup_pct":
            sens_sup,

        "Mejor_sensibilidad_individual_pct":
            mejor_sens_ind,
        "Ganancia_OR_vs_mejor_individual_pp":
            ganancia_or,
        "Ganancia_OR_bootstrap_IC95_inf_pp":
            gain_inf,
        "Ganancia_OR_bootstrap_IC95_sup_pp":
            gain_sup,
        "Ganancia_OR_bootstrap_mediana_pp":
            gain_mediana,
        "Ganancia_IC95_incluye_0":
            (
                "Sí"
                if (
                    not pd.isna(gain_inf)
                    and gain_inf <= 0 <= gain_sup
                )
                else (
                    "No"
                    if not pd.isna(gain_inf)
                    else "NA"
                )
            ),

        "N_B_LB_ClinVar":
            n_b_lb,
        "N_B_LB_falsos_positivos_por_ge1":
            n_b_lb_fp_or,
        "Falsos_positivos_OR_ge1_pct":
            fp_or,
        "Falsos_positivos_OR_IC95_inf_pct":
            fp_inf,
        "Falsos_positivos_OR_IC95_sup_pct":
            fp_sup,
        "Especificidad_OR_ge1_pct":
            especificidad_or,

        #Añado una columna de sensibilidad por herramienta
        "Sens_FATHMM-MKL_pct":
            sens_individual_salida(nombre_exacto="FATHH"),
        "Sens_UTRAnnotator_pct":
            sens_individual_salida(nombre_exacto="UTRannotator"),
        "Sens_GREEN-DB_pct":
            sens_individual_salida(nombre_exacto="GREEN-DB(3)"),
        "Sens_CADD_pct":
            sens_individual_salida(prefijo="CADD("),
        "Sens_GERP++_pct":
            sens_individual_salida(nombre_exacto="GERP(2)"),
        "Sens_SpliceAI_pct":
            sens_individual_salida(prefijo="Clasificación SpliceAI("),
        "Sens_dbscSNV-ADA_pct":
            sens_individual_salida(nombre_exacto="ADA_dbscSNV(0.6)"),
        "Sens_dbscSNV-RF_pct":
            sens_individual_salida(nombre_exacto="RF_dbscSNV(0.6)"),
        "Sens_PDIVAS_pct":
            sens_individual_salida(nombre_exacto="PDIVAS(0.082)"),
    }


#Repito el análisis para todos los grupos y umbrales

def analizar_regulacion(
    df,
    n_bootstrap,
    rng,
):
    filas = []

    for consecuencia in TARGET_CONSEQUENCES:
        for tipo in ["SNV", "INDEL"]:
            for umbral_cadd, cadd_col in CADD_CONFIG.items():

                herramientas = herramientas_regulacion(
                    consecuencia,
                    tipo,
                    cadd_col,
                )

                filas.append(
                    analizar_estrato(
                        df=df,
                        consecuencia=consecuencia,
                        tipo=tipo,
                        herramientas=herramientas,
                        umbral_cadd=umbral_cadd,
                        umbral_spliceai=None,
                        n_bootstrap=n_bootstrap,
                        rng=rng,
                    )
                )

    return pd.DataFrame(filas)


def analizar_splicing(
    df,
    n_bootstrap,
    rng,
):
    filas = []

    for consecuencia in TARGET_CONSEQUENCES:
        for tipo in ["SNV", "INDEL"]:
            for umbral_cadd, cadd_col in CADD_CONFIG.items():
                for (
                    umbral_spliceai,
                    spliceai_col,
                ) in SPLICEAI_CONFIG.items():

                    herramientas = herramientas_splicing(
                        consecuencia,
                        tipo,
                        cadd_col,
                        spliceai_col,
                    )

                    filas.append(
                        analizar_estrato(
                            df=df,
                            consecuencia=consecuencia,
                            tipo=tipo,
                            herramientas=herramientas,
                            umbral_cadd=umbral_cadd,
                            umbral_spliceai=umbral_spliceai,
                            n_bootstrap=n_bootstrap,
                            rng=rng,
                        )
                    )

    return pd.DataFrame(filas)


#Doy formato a las tablas del TFM

def numero_es(valor, decimales=2):
    if pd.isna(valor):
        return "NA"

    num = float(valor)

    texto = f"{num:,.{decimales}f}"
    return (
        texto
        .replace(",", "X")
        .replace(".", ",")
        .replace("X", ".")
    )


def entero_es(valor):
    if pd.isna(valor):
        return "NA"
    return f"{int(valor):,}".replace(",", ".")


def pct_n(porcentaje, n):
    if pd.isna(porcentaje):
        return "NA"

    return (
        f"{numero_es(porcentaje)} %\n"
        f"N = {entero_es(n)}"
    )


def pct_fraccion_ic(
    porcentaje,
    numerador,
    denominador,
    inf,
    sup,
):
    if pd.isna(porcentaje):
        return "NA"

    return (
        f"{numero_es(porcentaje)} %\n"
        f"{entero_es(numerador)}/{entero_es(denominador)}\n"
        f"IC95 % {numero_es(inf)}–{numero_es(sup)} %"
    )


def ganancia_ic(
    ganancia,
    inf,
    sup,
):
    if pd.isna(ganancia):
        return "NA"

    signo = "+" if float(ganancia) > 0 else ""

    return (
        f"{signo}{numero_es(ganancia)} pp\n"
        f"IC95 % bootstrap "
        f"{numero_es(inf)}–{numero_es(sup)} pp"
    )


def formatear_herramientas(texto):
    #Cambio los nombres para la tabla final
    if pd.isna(texto):
        return ""

    s = str(texto)

    reemplazos = {
        "FATHH": "FATHMM-MKL (0,5)",
        "GREEN-DB(3)": "GREEN-DB (3)",
        "GERP(2)": "GERP++ (2)",
        "UTRannotator": "UTRAnnotator",
        "CADD(10)": "CADD (10)",
        "CADD(15)": "CADD (15)",
        "CADD(20)": "CADD (20)",
        "Clasificación SpliceAI(0.2)": "SpliceAI (0,2)",
        "Clasificación SpliceAI(0.8)": "SpliceAI (0,8)",
        "ADA_dbscSNV(0.6)": "dbscSNV-ADA (0,6)",
        "RF_dbscSNV(0.6)": "dbscSNV-RF (0,6)",
        "PDIVAS(0.082)": "PDIVAS (0,082)",
    }

    for antiguo, nuevo in reemplazos.items():
        s = s.replace(antiguo, nuevo)

    return s


def tabla_principal(df):
    out = pd.DataFrame()

    out["Consecuencia"] = df["Consecuencia"]
    out["Tipo"] = df["Tipo"]
    out["Herramientas"] = (
        df["Herramientas_comparadas"]
        .apply(formatear_herramientas)
    )

    out["Cobertura conjunta"] = [
        pct_n(pct_, n)
        for pct_, n in zip(
            df["Cobertura_conjunta_pct"],
            df["N_evaluadas"],
        )
    ]

    out["Concordancia herramientas"] = [
        pct_n(pct_, n)
        for pct_, n in zip(
            df["Concordancia_herramientas_pct"],
            df["N_consenso_herramientas"],
        )
    ]

    out["Concordancia con ClinVar"] = [
        pct_n(pct_, n)
        for pct_, n in zip(
            df[
                "Concordancia_consenso_ClinVar_sobre_Nevaluadas_pct"
            ],
            df["N_consenso_correcto_ClinVar"],
        )
    ]

    out["Consenso P"] = [
        pct_n(pct_, n)
        for pct_, n in zip(
            df[
                "Consenso_patogenico_sobre_Nevaluadas_pct"
            ],
            df["N_consenso_patogenico"],
        )
    ]

    out["P→VUS"] = [
        pct_n(pct_, n)
        for pct_, n in zip(
            df["Consenso_P_ClinVar_VUS_pct"],
            df["Consenso_P_ClinVar_VUS_N"],
        )
    ]

    out["P→Conflicting"] = [
        pct_n(pct_, n)
        for pct_, n in zip(
            df["Consenso_P_ClinVar_Conflicting_pct"],
            df["Consenso_P_ClinVar_Conflicting_N"],
        )
    ]

    out["Sensibilidad ≥1 P"] = [
        pct_fraccion_ic(
            pct_,
            num,
            den,
            inf,
            sup,
        )
        for pct_, num, den, inf, sup in zip(
            df["Sensibilidad_OR_ge1_pct"],
            df["N_P_LP_detectadas_por_ge1"],
            df["N_P_LP_ClinVar"],
            df["Sensibilidad_OR_IC95_inf_pct"],
            df["Sensibilidad_OR_IC95_sup_pct"],
        )
    ]

    out["FP ≥1 P"] = [
        pct_fraccion_ic(
            pct_,
            num,
            den,
            inf,
            sup,
        )
        for pct_, num, den, inf, sup in zip(
            df["Falsos_positivos_OR_ge1_pct"],
            df["N_B_LB_falsos_positivos_por_ge1"],
            df["N_B_LB_ClinVar"],
            df["Falsos_positivos_OR_IC95_inf_pct"],
            df["Falsos_positivos_OR_IC95_sup_pct"],
        )
    ]

    out["Ganancia"] = [
        ganancia_ic(g, inf, sup)
        for g, inf, sup in zip(
            df["Ganancia_OR_vs_mejor_individual_pp"],
            df["Ganancia_OR_bootstrap_IC95_inf_pp"],
            df["Ganancia_OR_bootstrap_IC95_sup_pp"],
        )
    ]

    return out





#Programa principal

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Análisis comparativo final del TFM con IC95% Wilson "
            "y bootstrap pareado de la ganancia."
        )
    )

    parser.add_argument(
        "tabla_fusionada",
        help="Ruta a tabla_fusionada.tsv",
    )

    parser.add_argument(
        "--outdir",
        default="resultados_bootstrap",
        help=(
            "Directorio de salida "
            "(por defecto: resultados_bootstrap)."
        ),
    )

    parser.add_argument(
        "--n-bootstrap",
        type=int,
        default=10000,
        help=(
            "Número de remuestreos bootstrap "
            "(por defecto: 10000)."
        ),
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=20260927,
        help="Semilla aleatoria para reproducibilidad.",
    )

    args = parser.parse_args()

    entrada = Path(args.tabla_fusionada)
    outdir = Path(args.outdir)

    if not entrada.exists():
        sys.exit(
            f"ERROR: no existe el archivo: {entrada}"
        )

    if args.n_bootstrap <= 0:
        sys.exit(
            "ERROR: --n-bootstrap debe ser mayor que 0."
        )

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )


    df = pd.read_csv(
        entrada,
        sep="\t",
        dtype=str,
        low_memory=False,
    )



    df = preparar_datos(df)


    rng = np.random.default_rng(args.seed)


    #Guardo las tablas grandes con todos los umbrales
    #Si no se evalúa una herramienta aparece "-"
    out_reg_det = (
        outdir
        / "comparativa_regulacion_todos_umbral_con_IC.tsv"
    )
    out_spl_det = (
        outdir
        / "comparativa_splicing_todos_umbral_con_IC.tsv"
    )

    regulacion.to_csv(
        out_reg_det,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    splicing.to_csv(
        out_spl_det,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    #Para el TFM uso CADD 10 y SpliceAI 0,2
    reg_main = regulacion.loc[
        pd.to_numeric(
            regulacion["Umbral_CADD"],
            errors="coerce",
        ) == 10
    ].copy()

    spl_cadd = pd.to_numeric(
        splicing["Umbral_CADD"],
        errors="coerce",
    )
    spl_spliceai = pd.to_numeric(
        splicing["Umbral_SpliceAI"],
        errors="coerce",
    )

    spl_main = splicing.loc[
        (spl_cadd == 10)
        & (spl_spliceai == 0.2)
    ].copy()

    tabla_reg = tabla_principal(reg_main)
    tabla_spl = tabla_principal(spl_main)

    out_reg_main = (
        outdir
        / "tabla_TFM_regulacion_CADD10_con_IC.tsv"
    )
    out_spl_main = (
        outdir
        / "tabla_TFM_splicing_CADD10_SpliceAI02_con_IC.tsv"
    )

    tabla_reg.to_csv(
        out_reg_main,
        sep="\t",
        index=False,
    )

    tabla_spl.to_csv(
        out_spl_main,
        sep="\t",
        index=False,
    )


    print("\nArchivos generados:")
    for ruta in [
        out_reg_det,
        out_spl_det,
        out_reg_main,
        out_spl_main,
    ]:
        print(f"  {ruta}")





if __name__ == "__main__":
    main()