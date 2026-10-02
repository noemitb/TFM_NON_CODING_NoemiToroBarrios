#!/usr/bin/env python3

import pandas as pd

import numpy as np

import os

from statsmodels.stats.proportion import proportion_confint

# Noemí Toro Barrios | TFM: Revisión de herramientas y recursos para la anotación y priorización de variantes en regiones no codificantes.

#Script para calcular los parámetros estadísticos del análisis individual. 

# Archivo que contiene todas las herramientas ya fusionadas

INPUT_FILE = "tabla_fusionada.tsv"


# Nombre de los archivos donde guardaré los resultados

OUTPUT_DIR = "."

OUTPUT_FILE = "rendimiento_herramientas_grande.tsv"

OUTPUT_FILE_COMPACTO = "tabla_rendimiento_herramientas_nueva.tsv"


# Número de réplicas que voy a utilizar para calcular el intervalo de confianza del MCC mediante bootstrap

N_BOOTSTRAP = 10000


# Semilla para que el resultado del bootstrap sea reproducible si vuelvo a ejecutar el script

RANDOM_SEED = 19


# Columnas que quiero evaluar

HERRAMIENTAS = ["FATHH","UTRannotator","GREEN-DB(3)","CADD(10)","CADD(15)","CADD(20)","GERP(2)","Clasificación SpliceAI(0.2)","Clasificación SpliceAI(0.8)","ADA_dbscSNV(0.6)","RF_dbscSNV(0.6)","PDIVAS(0.082)"]


# Nombres que quiero utilizar después en la tabla final

NOMBRES_HERRAMIENTAS = {"FATHH": "FATHMM-MKL (0,5)","UTRannotator": "UTRAnnotator","GREEN-DB(3)": "GREEN-DB (3)","CADD(10)": "CADD (10)","CADD(15)": "CADD (15)","CADD(20)": "CADD (20)","GERP(2)": "GERP++ (2)","Clasificación SpliceAI(0.2)": "SpliceAI (0,2)","Clasificación SpliceAI(0.8)": "SpliceAI (0,8)","ADA_dbscSNV(0.6)": "dbscSNV-ADA (0,6)","RF_dbscSNV(0.6)": "dbscSNV-RF (0,6)","PDIVAS(0.082)": "PDIVAS (0,082)"}


# Categorías de ClinVar que voy a utilizar como referencia para calcular sensibilidad, especificidad, VPP y MCC
#En posteriores versiones, arreglar el PBenigna.
CLINVAR_MAP = {"PBenigna/Probablemente benigna": "Benigna","Benigna/Probablemente benigna": "Benigna","Patogénica/Probablemente patogénica": "Patogénica"}


# Algunas herramientas tienen las predicciones escritas de una forma diferente, así que las dejo todas iguales(arreglar misma salida para todos en un futuro)

MAP_HERRAMIENTAS = {"patogénico": "Patogénica","patogénica": "Patogénica","probablemente patogénico": "Patogénica","probablemente patogénica": "Patogénica","benigno": "Benigna","benigna": "Benigna","probablemente benigno": "Benigna","probablemente benigna": "Benigna","vus": np.nan}


# Valores que significan que la herramienta no ha podido anotar la variante

NA_VALUES = ["No anotada","No anotado","No concluyente"]


# Consecuencias codificantes que NO INCLUYO

CONSECUENCIAS_CODIFICANTES = ["missense_variant","synonymous_variant","nonsense","frameshift_variant","inframe_deletion","inframe_insertion","inframe_indel","stop_gained","stop_lost","start_lost","initiator_codon_variant","protein_altering_variant","coding_sequence_variant"]


# Calculo el porcentaje y el intervalo de confianza de Wilson

def calcular_ic_wilson(n_exitos, n_total):


    # Si no tengo variantes para calcular la métrica, dejo el resultado como NaN

    if n_total == 0:

        return np.nan, np.nan, np.nan


    porcentaje = (n_exitos / n_total) * 100


    lim_inf, lim_sup = proportion_confint(count=n_exitos,nobs=n_total,alpha=0.05,method="wilson" )


    return (porcentaje,lim_inf * 100,lim_sup * 100)


# Compruebo si una variante tiene alguna consecuencia codificante

def tiene_consecuencia_codificante(consecuencia):


    if pd.isna(consecuencia):

        return False


    if consecuencia == "Sin_Dato":

        return False


    consecuencias = str(consecuencia).split("|")


    return any(

        cons.strip() in CONSECUENCIAS_CODIFICANTES

        for cons in consecuencias

    )


# Calculo el coeficiente de correlación de Matthews utilizando los cuatro valores de la matriz de confusión

def calcular_mcc(VP, VN, FP, FN):


    denominador = np.sqrt((VP + FP)* (VP + FN)* (VN + FP)* (VN + FN))


    # Si el denominador es cero, el MCC no se puede calcular

    if denominador == 0:

        return np.nan


    return ((VP * VN) - (FP * FN)) / denominador


# Calculo un intervalo de confianza del MCC mediante bootstrap remuestreando las variantes con reemplazo 10000 veces

def bootstrap_mcc(VP, VN, FP, FN, rng):


    total = (VP + VN + FP + FN)


    # Si no tengo variantes,no puedo realizar el bootstrap

    if total == 0:

        return np.nan, np.nan


    # Cada variante evaluada pertenece a una de las cuatro categorías de la matriz de confusión: VP, VN, FP o FN

    probabilidades = np.array([VP / total,VN / total, FP / total,FN / total])


    # Para el MCC solamente necesito conocer cuántas variantes caen en cada una de estas cuatro categorías
    # Este remuestreo es equivalente a seleccionar variantes con reemplazo una a una, pero es mucho más rápido

    replicas = rng.multinomial(total,probabilidades,size=N_BOOTSTRAP )


    VP_boot = replicas[:, 0]

    VN_boot = replicas[:, 1]

    FP_boot = replicas[:, 2]

    FN_boot = replicas[:, 3]


    denominador = np.sqrt((VP_boot + FP_boot)* (VP_boot + FN_boot)* (VN_boot + FP_boot)* (VN_boot + FN_boot) )


    # Algunas réplicas, especialmente si el número de variantes es pequeño, podrían tener un denominador igual a cero

    mcc_boot = np.full( N_BOOTSTRAP, np.nan )


    validas = ( denominador > 0)


    mcc_boot[validas] = ((VP_boot[validas] * VN_boot[validas] - FP_boot[validas] * FN_boot[validas]) / denominador[validas])


    # Utilizo el método percentil:

    # percentil 2,5 y percentil 97,5 del bootstrap

    lim_inf = np.nanpercentile(mcc_boot, 2.5)


    lim_sup = np.nanpercentile( mcc_boot,97.5)


    return ( lim_inf,lim_sup)


def main():


    # Cargo la tabla fusionada

    print(f"Cargando archivo: {INPUT_FILE}...")


    df = pd.read_csv(INPUT_FILE,sep="\t",low_memory=False)


    # Compruebo que todas las herramientas que quiero analizar estén presentes en la tabla

    for herramienta in HERRAMIENTAS:


        if herramienta not in df.columns:

            print(f"Error: no encuentro la columna "f"'{herramienta}' en la tabla." )

            return


    # Quito las consecuencias codificantes y también las variantes sin consecuencia conocida

    


    df = df[(df["Consecuencia"].notna()) &  (~df["Consecuencia"].apply(tiene_consecuencia_codificante )) &(df["Consecuencia"] != "Sin_Dato")].copy()


    # Convierto las variantes no anotadas en NaN
        

    df[HERRAMIENTAS] = df[HERRAMIENTAS].replace(NA_VALUES,np.nan)


    # Homogeneizo las categorías que utilizan algunas herramientas

    df[HERRAMIENTAS] = df[HERRAMIENTAS].replace(MAP_HERRAMIENTAS)


    # Convierto las categorías de ClinVar en las dos clases

   
    df["ClinVar"] = df["Clasificación Clinvar"].map(CLINVAR_MAP)


    # Me quedo solamente con variantes de ClinVar patogénicas o benignas, VUS, Conflicting y Otras no se utilizan para estas métricas

    df_eval = df[df["ClinVar"].notna()].copy()


    resultados = []


    # Creo el generador aleatorio que utilizaré para que el bootstrap sea reproducible (para recordarlo: al crear la semilla garantizo que se cree la misma aletoriedad en futura reproducciones del script y no me vuelva loca porque me sale resultados distintos)


    rng = np.random.default_rng(RANDOM_SEED)


    # Calculo las métricas herramienta por herramienta

    for herramienta in HERRAMIENTAS:


        # Cada herramienta tiene una cobertura diferente, así que utilizo únicamente las variantes para las que ha generado una predicción válida

        sub = df_eval[df_eval[herramienta].isin(["Patogénica", "Benigna"])].copy()


        # ClinVar me dice cuál es la clasificación de referencia

        y_true = (sub["ClinVar"] == "Patogénica")


        # La herramienta me dice cuál es su predicción

        y_pred = (sub[herramienta] == "Patogénica")


        # Verdaderos positivos:

        # ClinVar patogénica y herramienta patogénica

        VP = (y_true & y_pred).sum()


        # Verdaderos negativos:

        # ClinVar benigna y herramienta benigna

        VN = (~y_true & ~y_pred).sum()


        # Falsos positivos:

        # ClinVar benigna pero herramienta patogénica

        FP = (~y_true & y_pred).sum()


        # Falsos negativos:

        # ClinVar patogénica pero herramienta benigna

        FN = (y_true & ~y_pred).sum()


        # Número total de variantes evaluadas por la herramienta

        total = (VP + VN + FP + FN)


        # Total de variantes patogénicas según ClinVar

        total_patos = (VP + FN)


        # Total de variantes benignas según ClinVar

        total_benignas = (VN + FP)


        # Total de variantes que la herramienta predice como patogénicas

        total_patos_herramienta = (VP + FP)


        # Prevalencia de variantes P/LP:

        
        if total > 0:

            prevalencia = (total_patos / total) * 100

        else:

            prevalencia = np.nan


        # Sensibilidad (de todas las patogénicas de ClinVar,cuántas detecta correctamente la herramienta)

        sensibilidad, sen_inf, sen_sup = calcular_ic_wilson(VP,total_patos)


        # Especificidad(de todas las benignas de ClinVar,cuántas identifica correctamente como benignas)

        especificidad, esp_inf, esp_sup = calcular_ic_wilson(VN,total_benignas)


        # VPP(de todas las que la herramienta llama patogénicas, cuántas son realmente patogénicas según ClinVar)

        vpp, vpp_inf, vpp_sup = calcular_ic_wilson(VP,total_patos_herramienta)


        # MCC (medida para tener en cuenta el desequilibrio de clases, -1 completamente inverso a lo esperado, 0 al azar, 1 perfecto)

        mcc = calcular_mcc(VP,VN,FP,FN)


        # Calculo el intervalo de confianza del MCC mediante 10000 réplicas bootstrap

        mcc_inf, mcc_sup = bootstrap_mcc(VP,VN,FP,FN,rng)

        # Guardo los resultados de esta herramienta

        resultados.append({
            "Herramienta":

                NOMBRES_HERRAMIENTAS[herramienta],


            "N_evaluadas":

                total,


            "N_P/LP":

                total_patos,


            "Prevalencia_P/LP_%":

                round(prevalencia, 2),


            "Sensibilidad_%":

                round(sensibilidad, 2),


            "Sensibilidad_IC95":

                f"[{sen_inf:.2f} - {sen_sup:.2f}]",


            "Especificidad_%":

                round(especificidad, 2),


            "Especificidad_IC95":

                f"[{esp_inf:.2f} - {esp_sup:.2f}]",


            "VPP_%":

                round(vpp, 2),


            "VPP_IC95":

                f"[{vpp_inf:.2f} - {vpp_sup:.2f}]",


            "MCC":

                round(mcc, 3),


            "MCC_IC95_bootstrap":

                f"[{mcc_inf:.3f} - {mcc_sup:.3f}]",


            "VP":

                VP,


            "VN":

                VN,


            "FP":

                FP,


            "FN":

                FN

        })


    # Paso todos los resultados a una tabla

    resultados_df = pd.DataFrame(resultados )


    # Como no cabe todo en el TFM, guardo también una tabla más pequeña

    resultados_compactos = resultados_df[[

        "Herramienta",

        "N_evaluadas",

        "N_P/LP",

        "Prevalencia_P/LP_%",

        "Sensibilidad_%",

        "Sensibilidad_IC95",

        "Especificidad_%",

        "Especificidad_IC95",

        "VPP_%",

        "VPP_IC95",

        "MCC",

        "MCC_IC95_bootstrap"

    ]].copy()


    # Guardo los archivos finales

    os.makedirs(

        OUTPUT_DIR,

        exist_ok=True

    )


    ruta_salida = os.path.join(

        OUTPUT_DIR,

        OUTPUT_FILE

    )


    ruta_salida_compacta = os.path.join(

        OUTPUT_DIR,

        OUTPUT_FILE_COMPACTO

    )


    resultados_df.to_csv(

        ruta_salida,

        sep="\t",

        index=False

    )


    resultados_compactos.to_csv(

        ruta_salida_compacta,

        sep="\t",

        index=False

    )

if __name__ == "__main__":

    main()
