#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import sys
import numpy as np

# Noemí Toro Barrios | TFM - Revisión de herramientas y recursos para la anotación y priorización de variantes en regiones no codificantes
#Generar HeatMaps


def procesar_tabla(input_file):

    # Cargo la tabla
    print(f"Cargando archivo: {input_file}...")
    df = pd.read_csv(input_file, sep="\t")


    # Selecciono las herramientas de splicing
    tools = [
        "Clasificación SpliceAI(0.2)",
        "Clasificación SpliceAI(0.8)",
        "ADA_dbscSNV(0.6)",
        "RF_dbscSNV(0.6)",
        "PDIVAS(0.082)",
        "GERP(2)",
        "CADD(10)",
        "CADD(15)",
        "CADD(20)"
    ]


    # Compruebo que todas las herramientas estén en la tabla
    for tool in tools:
        if tool not in df.columns:
            print(f"Error: no encuentro la columna '{tool}'")
            sys.exit(1)


    # Paso los valores sin anotación a NaN
    # para que no entren en el denominador
    NA_VALUES = [
        "No anotada",
        "No anotado",
        "No concluyente"
    ]

    df[tools] = df[tools].replace(
        NA_VALUES,
        np.nan
    )


    # Quito las consecuencias codificantes porque quiero estudiar
    # únicamente variantes de regiones no codificantes
    consecuencias_codificantes = [
        "missense_variant",
        "synonymous_variant",
        "nonsense",
        "frameshift_variant",
        "inframe_deletion",
        "inframe_insertion",
        "inframe_indel",
        "stop_gained",
        "stop_lost",
        "start_lost",
        "initiator_codon_variant",
        "protein_altering_variant",
        "coding_sequence_variant"
    ]


    print("Filtrando variantes...")


    def tiene_consecuencia_codificante(consecuencia):

        if pd.isna(consecuencia) or consecuencia == "Sin_Dato":
            return False

        consecuencias = str(consecuencia).split("|")

        return any(
            cons.strip() in consecuencias_codificantes
            for cons in consecuencias
        )


    df_no_codificante = df[
        (~df["Consecuencia"].apply(tiene_consecuencia_codificante)) &
        (df["Consecuencia"] != "Sin_Dato") &
        (df["Clasificación Clinvar"] != "Otras")
    ].copy()


    # Separo las consecuencias múltiples después de eliminar
    # las variantes que presentan alguna consecuencia codificante
    df_no_codificante["Consecuencia"] = (
        df_no_codificante["Consecuencia"]
        .astype(str)
        .str.split("|")
    )

    df_no_codificante = df_no_codificante.explode(
        "Consecuencia"
    )

    df_no_codificante["Consecuencia"] = (
        df_no_codificante["Consecuencia"]
        .str.strip()
    )


    # Defino las distintas formas en las que puede aparecer
    # una predicción patogénica
    terminos_patogenicos = [
        "patogénica",
        "patogénico",
        "probablemente patogénico",
        "probablemente patogénica"
    ]


    # Esta tabla contendrá los porcentajes usados en el heatmap
    results = []

    # Esta tabla contendrá únicamente el número de variantes
    # con una predicción válida para material suplementario
    results_n = []


    for (cons, clin), group in df_no_codificante.groupby(
        ["Consecuencia", "Clasificación Clinvar"]
    ):

        # Fila para el heatmap
        row = {
            "Consecuencia": cons,
            "ClinVar": clin
        }

        # Fila para la tabla suplementaria
        row_n = {
            "Consecuencia": cons,
            "ClinVar": clin
        }


        # Cada herramienta puede haber anotado un número
        # diferente de variantes
        for tool in tools:

            # Me quedo solo con las variantes con
            # predicción válida
            valores_tool = group[tool].dropna()

            # N = número de variantes anotadas
            N = len(valores_tool)

            # Guardo N en la tabla suplementaria
            if N == 0:

                perc = np.nan
                row_n[tool] = "—"

            else:

                valores_tool = (
                    valores_tool
                    .astype(str)
                    .str.strip()
                    .str.lower()
                )

                # Identifico las variantes que superan el umbral
                es_patogenica = valores_tool.isin(
                    terminos_patogenicos
                )

                # n = número que supera el umbral
                n = es_patogenica.sum()

                # Porcentaje para el heatmap
                perc = (n / N) * 100

                # Tabla suplementaria: solo n/N
                row_n[tool] = f"{n}/{N}"


            if N == 0:

                # Si no hay variantes anotadas,
                # dejo la celda vacía en el heatmap
                perc = np.nan

            else:

                # Normalizo los textos
                valores_tool = (
                    valores_tool
                    .astype(str)
                    .str.strip()
                    .str.lower()
                )

                # Calculo el porcentaje de variantes anotadas
                # que supera el umbral de patogenicidad
                perc = (
                    valores_tool
                    .isin(terminos_patogenicos)
                    .mean()
                    * 100
                )


            row[tool] = perc


        results.append(row)
        results_n.append(row_n)


    # Si no quedan datos después del filtrado
    if not results:
        print("Error: No hay datos suficientes.")
        sys.exit(1)


    # Tabla de porcentajes
    df_results = pd.DataFrame(results)

    df_results = df_results.sort_values(
        by=["Consecuencia", "ClinVar"]
    ).reset_index(drop=True)


    # Tabla con el número de variantes anotadas
    df_n = pd.DataFrame(results_n)

    df_n = df_n.sort_values(
        by=["Consecuencia", "ClinVar"]
    ).reset_index(drop=True)


    # Me quedo solo con las herramientas
    # para crear el heatmap
    heatmap_df = df_results[tools]


    print(
        "Generando mapa de calor horizontal "
        "de herramientas de splicing..."
    )


    # Ajusto el tamaño del gráfico al número
    # de herramientas y filas
    ancho = max(
        12,
        len(tools) * 1.5
    )

    alto = max(
        10,
        len(df_results) * 0.45
    )


    fig, ax = plt.subplots(
        figsize=(ancho, alto)
    )


    # Creo el heatmap
    sns.heatmap(
        heatmap_df,
        annot=True,
        cmap="YlOrRd",
        fmt=".1f",
        cbar_kws={
            "label": "% Variantes que superan el umbral"
        },
        ax=ax,
        linewidths=0.5,
        linecolor="lightgray",
        vmin=0,
        vmax=100
    )


    # Coloco las categorías de ClinVar en el eje Y
    ax.set_yticks(
        [i + 0.5 for i in range(len(df_results))]
    )


    etiquetas_clinvar = (
        df_results["ClinVar"]
        .replace({
            "PBenigna/Probablemente benigna":
                "Benigna/Prob. benigna",

            "Patogénica/Probablemente patogénica":
                "Patogénica/Prob. patogénica"
        })
    )


    ax.set_yticklabels(
        etiquetas_clinvar,
        rotation=0,
        fontsize=10
    )

    ax.set_ylabel("")

    ax.tick_params(
        axis="y",
        length=0,
        pad=10
    )


    # Colores para distinguir visualmente
    # las categorías clínicas de ClinVar
    colores_clinvar = {
        "Patogénica/Probablemente patogénica": "#ffcccc",
        "PBenigna/Probablemente benigna": "#ccffcc",
        "VUS": "#fff2cc",
        "Conflicting": "#f2f2f2"
    }


    # Este límite controla hasta dónde se extiende hacia
    # la izquierda el fondo de color de ClinVar
    limite_izquierdo = -0.45


    # Nombres preparados manualmente para que las consecuencias
    # largas no se monten unas encima de otras
    nombres_consecuencias = {

        "3_prime_UTR_variant":
            "3 prime UTR\nvariant",

        "5_prime_UTR_variant":
            "5 prime UTR\nvariant",

        "genic_downstream_transcript_variant":
            "Genic downstream\ntranscript variant",

        "genic_upstream_transcript_variant":
            "Genic upstream\ntranscript variant",

        "intron_variant":
            "Intron\nvariant",

        "no_sequence_alteration":
            "No sequence\nalteration",

        "non_coding_transcript_variant":
            "Non-coding transcript\nvariant",

        "non-coding-transcript_variant":
            "Non-coding transcript\nvariant",

        "splice_acceptor_variant":
            "Splice acceptor\nvariant",

        "splice_donor_variant":
            "Splice donor\nvariant"
    }


    consecuencias = df_results[
        "Consecuencia"
    ].unique()


    # Coloco a la izquierda el nombre de cada consecuencia
    # y separo los bloques con líneas negras
    for cons in consecuencias:

        indices = df_results.index[
            df_results["Consecuencia"] == cons
        ].tolist()

        inicio = indices[0]
        fin = indices[-1] + 1

        centro = (
            inicio + fin
        ) / 2.0


        if cons in nombres_consecuencias:

            cons_limpia = nombres_consecuencias[cons]

        else:

            cons_limpia = (
                cons
                .replace("_", " ")
                .capitalize()
            )


        # Coloco la consecuencia justo a la izquierda
        # del bloque de colores de ClinVar
        ax.text(
            limite_izquierdo - 0.03,
            centro,
            cons_limpia,
            ha="right",
            va="center",
            fontweight="bold",
            fontsize=11,
            color="#333333",
            linespacing=1.0,
            multialignment="right",
            transform=ax.get_yaxis_transform()
        )


        # Separo cada consecuencia con una línea negra
        if inicio > 0:

            ax.axhline(
                inicio,
                xmin=limite_izquierdo,
                xmax=1,
                color="black",
                linewidth=2.5,
                clip_on=False
            )


    for i in range(len(df_results)):

        clin_val = df_results.loc[
            i,
            "ClinVar"
        ]

        color_fondo = colores_clinvar.get(
            clin_val,
            "#ffffff"
        )


        ax.axhspan(
            i,
            i + 1,
            xmin=limite_izquierdo,
            xmax=0,
            facecolor=color_fondo,
            alpha=1.0,
            clip_on=False,
            zorder=-1
        )


    nombres_tools = {

        "Clasificación SpliceAI(0.2)":
            "SpliceAI (0,2)",

        "Clasificación SpliceAI(0.8)":
            "SpliceAI (0,8)",

        "ADA_dbscSNV(0.6)":
            "dbscSNV-ADA (0,6)",

        "RF_dbscSNV(0.6)":
            "dbscSNV-RF (0,6)",

        "PDIVAS(0.082)":
            "PDIVAS (0,082)",

        "GERP(2)":
            "GERP++ (2)",

        "CADD(10)":
            "CADD (10)",

        "CADD(15)":
            "CADD (15)",

        "CADD(20)":
            "CADD (20)"
    }


    ax.set_xticklabels(
        [
            nombres_tools.get(tool, tool)
            for tool in tools
        ],
        rotation=45,
        ha="right",
        fontsize=10
    )


    # Cambio los nombres de las herramientas
    # en la tabla suplementaria
    df_n = df_n.rename(
        columns=nombres_tools
    )


    # Cambio las categorías de ClinVar
    # para que sean más fáciles de leer
    df_n["ClinVar"] = df_n["ClinVar"].replace({

        "PBenigna/Probablemente benigna":
            "Benigna/Prob. benigna",

        "Patogénica/Probablemente patogénica":
            "Patogénica/Prob. patogénica"
    })


    plt.title(
        "Herramientas de splicing",
        fontsize=16,
        pad=20
    )


    # Dejo margen suficiente para las consecuencias
    # y categorías de ClinVar
    plt.subplots_adjust(
        left=0.42
    )


    # Guardo la figura
    output_img = "heatmap_horizontal_splicing.png"

    plt.savefig(
        output_img,
        bbox_inches="tight",
        dpi=300
    )


    # Guardo la tabla destinada al material suplementario
    df_n.to_csv(
        "tabla_suplementaria_N_splicing.tsv",
        sep="\t",
        index=False
    )


    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    args = parser.parse_args()

    procesar_tabla(args.input)