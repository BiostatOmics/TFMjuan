# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 13:31:32 2024

Compositional Analysis

@author: Juan Andrés Tejedor Serrano
"""

# ============================================
# CARGA DE MÓDULOS
# ============================================
import argparse
import pandas as pd
import anndata as ad
import scanpy as sc
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz
from sccoda.util import comp_ana as mod
import matplotlib.pyplot as plt
import arviz as az

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.curdir, '02-scripts')))
from my_utils import log_message
# ============================================
# PARSING DE CLI A PYTHON
# ============================================

# Creo el objeto parseador
parser = argparse.ArgumentParser(description="Plots")

# Establezco argumentos a parsear
parser.add_argument("input_amp", type=str)
parser.add_argument("results_path", type=str)

# Extraigo los argumentos del CLI a un objeto de Python
args = parser.parse_args()

# Recupero argumentos como variables
input_amp = args.input_amp
results_path = args.results_path


# ============================================
# CONFIGURACIÓN INICIAL Y CARGA DE DATOS
# ============================================

# Establezco ajustes de Scanpy
sc.settings.verbosity = 4
sc.settings.n_jobs = 8

# Cargo datos y veo sus atributos
log_message(f"Cargando amp desde {input_amp}...")
amp = ad.read_h5ad(input_amp)
log_message(f"Datos cargados desde {input_amp}:\n")
print(amp)

# Filtro y elimino las células Unknown
amp = amp[amp.obs["annotations"] != "Unknown"].copy()
# ============================================
# DEFINICIÓN DE FUNCIONES
# ============================================


def ensure_directory_exists(directory):
    # Crea el directorio especificado en el caso de que no exista
    if not os.path.exists(directory):
        os.makedirs(directory)


def convert_format(adata, region):
    # Filtrar por la región cerebral si se especifica
    if region != "all":
        adata = adata[adata.obs['brain_region'] == region]
    # ============================================
    # CONVERSIÓN DE FORMATO SCANPY A SCCODA
    # ============================================
    log_message("Convirtiendo formato...")
    # Obtengo un df sólo con los participantes (únicos) y su condición
    covariates_df = adata.obs[["participant_id", "case_control",
                               "path_braak_lb", "path_braak_nft",
                               "hoehn_yahr_stage"]].drop_duplicates()
    # Establezco el participante como índice
    covariates_df.set_index("participant_id", inplace=True)
    # Establezco el participante como índice
    coda = dat.from_scanpy(
        adata,
        cell_type_identifier="annotations",
        sample_identifier="participant_id",
        covariate_df=covariates_df
    )
    # Elimino sujetos "other"
    coda_case_control = coda[coda.obs["case_control"].isin(["Control", "Case"])]
    # Convierto a categoría
    coda_case_control.obs.path_braak_lb = coda_case_control.obs.path_braak_lb.astype("category")
    coda_case_control.obs.path_braak_nft = coda_case_control.obs.path_braak_nft.astype("category")
    coda_case_control.obs.hoehn_yahr_stage = coda_case_control.obs.hoehn_yahr_stage.astype("category")

    return coda_case_control


def generate_exploratory_plots(coda_case_control, results_path, region):
    # ============================================
    # GENERACIÓN DE PLOTS EXPLORATIVOS PRE-ANÁLISIS
    # ============================================
    
    features = ["case_control", "path_braak_lb",
                "path_braak_nft", "hoehn_yahr_stage"]
    
    for feature in features:
        log_message(f"Generando barplots para {feature}")
        # Stacked barplots
        plot_type = "barplot"
        base_file_name = f"{results_path}amppd_{plot_type}_{region}"
        viz.stacked_barplot(coda_case_control, feature_name=feature)
        plt.savefig(f"{base_file_name}_{feature}.png", bbox_inches="tight")
        plt.close()

        log_message(f"Generando boxplots para {feature}")
        # Grouped boxplots
        plot_type = "boxplot"
        base_file_name = f"{results_path}amppd_{plot_type}_{region}"
        viz.boxplots(
            coda_case_control,
            feature_name=feature,
            plot_facets=False,
            y_scale="relative",
            add_dots=False,
        )
        plt.savefig(f"{base_file_name}_{feature}.png", bbox_inches="tight")
        plt.close()

        log_message(f"Generando boxplots separados para {feature}")
        # Grouped boxplots per cell type
        plot_type = "boxplot_separated"
        base_file_name = f"{results_path}amppd_{plot_type}_{region}"
        viz.boxplots(
            coda_case_control,
            feature_name=feature,
            plot_facets=True,
            y_scale="log",
            add_dots=True,
            cmap="Reds",
        )
        plt.savefig(f"{base_file_name}_{feature}.png", bbox_inches="tight")
        plt.close()


def run_compositional_analysis(coda_case_control, results_path, region, fdr=0.05):
    log_message("Generando modelo composicional...")
    # Generación de modelo composicional (establezco "Control" como la categoría de referencia)
    model_case_control = mod.CompositionalAnalysis(coda_case_control,
                                                   formula="C(case_control, Treatment('Control'))",
                                                   reference_cell_type="automatic")

    # Ejecución de simulación para la inferencia de los parámetros vía MonteCarlo Sampling
    log_message("Ejecutando simulación para inferencia de parámetros...")
    results = model_case_control.sample_hmc()

    # Guardar resultados
    log_message("Guardando resultados del análisis composicional...")
    intercept_file = f"{results_path}intercept_case_control_{fdr}_{region}.tsv"
    effect_file = f"{results_path}effect_case_control_{fdr}_{region}.tsv"
    credible_effects_file = f"{results_path}credible_effects_case_control_{fdr}_{region}.tsv"
    
    results.intercept_df.to_csv(intercept_file, sep="\t", header=True, index=True)
    results.effect_df.to_csv(effect_file, sep="\t", header=True, index=True)
    results.credible_effects().to_csv(credible_effects_file, sep="\t", header=True, index=True)

    results.save(f"{results_path}scCODA_case_control_{region}.pkl")
    
    log_message(f"Resultados *EXTENDIDOS* del análisis composicional (FDR: {fdr}):\n")
    results.summary_extended()

    log_message(f"Resultados *SIGNIFICATIVOS* del análisis composicional (FDR: {fdr}):\n")
    print(results.credible_effects())
    
    return results


def generate_post_analysis_plots(results, results_path, region):
    log_message(f"Generando plots de inferencia bayesiana")
    plot_type = "bayesian_results"
    base_file_name = f"{results_path}amppd_{plot_type}_{region}"
    # Plots bayesianos
    az.plot_trace(
        results,
        divergences=False,
        var_names=["alpha", "beta"],
        coords={"cell_type": results.posterior.coords["cell_type_nb"]},
    )
    plt.savefig(f"{base_file_name}.png", bbox_inches="tight")
    plt.close()

    log_message(f"Generando plot de log2FC")
    plot_type = "log2FC"
    base_file_name = f"{results_path}amppd_{plot_type}_{region}"
    results.effect_df[["log2-fold change"]].reset_index(level=0, drop=True).plot.bar()
    plt.savefig(f"{base_file_name}.png", bbox_inches="tight")
    plt.close()


def compositional_analysis(amp, region, results_path):
    # Crear subcarpeta específica para la región
    region_results_path = f"{results_path}{region}/"
    
    ensure_directory_exists(region_results_path)

    # Paso 1: Conversión de formato
    amp_coda_case_control = convert_format(amp, region)
    log_message(f"Conversión finalizada:\n{amp_coda_case_control}\n{amp_coda_case_control.var}")

    # Paso 2: Generación de plots exploratorios pre-análisis
    generate_exploratory_plots(amp_coda_case_control, region_results_path, region)
    log_message("Plots exploratorios generados")

    # Paso 3: Análisis composicional casos vs controles
    results = run_compositional_analysis(amp_coda_case_control, region_results_path, region)
    log_message("Análisis composicional completado")

    # Paso 4: Generación de plots post-análisis
    generate_post_analysis_plots(results, region_results_path, region)
    log_message("Plots post-análisis generados")

    # Paso 5: Resultados adicionales con FDR = 20%
    fdr = 0.2
    results.set_fdr(est_fdr=fdr)
    log_message(f"Resultados *EXTENDIDOS* del análisis composicional (FDR: {fdr}):\n")
    results.summary_extended()

    log_message(f"Resultados *SIGNIFICATIVOS* del análisis composicional (FDR: {fdr}):\n")
    print(results.credible_effects())

    log_message("Guardando resultados del análisis composicional...")
    results.credible_effects().to_csv(f"{region_results_path}scCODA_credible_effects_case_control_{fdr}_{region}.tsv",
                                      sep="\t", header=True, index=True)


# Llamada a la función principal
compositional_analysis(amp=amp, region="all", results_path=results_path)
compositional_analysis(amp=amp, region="DMNX", results_path=results_path)
compositional_analysis(amp=amp, region="GPI", results_path=results_path)
compositional_analysis(amp=amp, region="PFC", results_path=results_path)
compositional_analysis(amp=amp, region="PMC", results_path=results_path)
compositional_analysis(amp=amp, region="PVC", results_path=results_path)
