# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 12:53:03 2024

Draft Plots for Exploration

@author: Juan Andrés Tejedor Serrano
"""

# ============================================
# CARGA DE MÓDULOS
# ============================================
import argparse
import numpy as np
import anndata as ad
import scanpy as sc
import symphonypy as sp
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

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
parser.add_argument("save_path", type=str)
# Extraigo los argumentos del CLI a un objeto de Python
args = parser.parse_args()

# Recupero argumentos como variables
input_amp = args.input_amp
save_path = args.save_path

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

# ============================================
# CONFIGURACIÓN DE PLOTS
# ============================================

sc.settings.autoshow = False  # No mostrar gráfico en consola. Permite realizar modificaciones a los gráficos obtenidos por scanpy para después guardarlos de forma explícita.
dpi_save = 150
sc.set_figure_params(dpi=80,
                     dpi_save=dpi_save,  # DPI de 300 para alta calidad
                     figsize=(14, 10),  # Tamaño de 14 pulgadas x 10 pulgadas
                     frameon=True,
                     color_map="magma")  # Mostrar marco (ejes)

xkcd_colors = list(matplotlib.colors.XKCD_COLORS.values())
plt.rcParams.update({
    'axes.titlesize': 24,     # Título del gráfico
    'axes.labelsize': 18,     # Etiquetas de los ejes
    'xtick.labelsize': 14,    # Texto en el eje X
    'ytick.labelsize': 14,    # Texto en el eje Y
    'legend.fontsize': 16,    # Tamaño de la fuente de la leyenda
    'font.size': 16,           # Tamaño de fuente general
})
size = 0.5


# ============================================
# GENERACIÓN DE PLOTS
# ============================================

# Reordeno aleatoriamente las células, evitando solapamiento de puntos.
np.random.seed(0)
random_indices = np.random.permutation(list(range(amp.shape[0])))
amp_random_order = amp[random_indices, :]

# Establezco el nombre base de los archivos
plot_type = "tsne_exploratory"
base_file_name = f"{save_path}amppd_{plot_type}"

categoric_to_plot = ["supercluster_term",
                     "cell_type",
                     "ROIGroupFine",
                     "leiden_1_00_corrected",
                     "leiden_0_25_corrected",
                     "leiden_0_10_corrected",
                     "set",
                     "case_control",
                     "brain_region",
                     "Phase",
                     ]

log_message("Generando tsne plots de variables categóricas...")
for color in categoric_to_plot:
    ncat = len(amp.obs[color].cat.categories)
    fig, ax = plt.subplots()
    if ncat <= 30:
        sc.pl.tsne(amp_random_order,
                   color=color,
                   ax=ax,
                   show=False,
                   s=size,
                   alpha=1,
                   palette=sns.husl_palette(ncat),
                   frameon=True,
                   )
    else:
        sc.pl.tsne(amp_random_order,
                   color=color,
                   ax=ax,
                   show=False,
                   s=size,
                   alpha=1,
                   palette=xkcd_colors,
                   frameon=True,
                   )
    if ncat <= 30:
        plt.legend(ncol=1, bbox_to_anchor=(1, 1), frameon=False)
    elif ncat > 30 and ncat <= 60:
        plt.legend(ncol=2, bbox_to_anchor=(1, 1), frameon=False)
    elif ncat > 60 and ncat <= 90:
        plt.legend(ncol=3, bbox_to_anchor=(1, 1), frameon=False)
    elif ncat > 90 and ncat <= 150:
        plt.legend(ncol=4, bbox_to_anchor=(1, 1), frameon=False)
    else:
        plt.legend(ncol=6, bbox_to_anchor=(1, 1), fontsize=6, frameon=False)
    file_name = f"{base_file_name}_{color}"
    plt.savefig(f"{file_name}.png", bbox_inches="tight")
    plt.close()

continous_to_plot = ["log1p_n_genes_by_counts",
                     "log1p_total_counts",
                     "doublet_score",
                     "pct_counts_in_top_20_genes",
                     "hoehn_yahr_stage",
                     "path_braak_lb",
                     "path_braak_nft",
                     "percent_apop",
                     "percent_dna_repair",
                     "percent_ieg",
                     "percent_mito",
                     "percent_mito_ribo",
                     "percent_oxphos",
                     "symphony_per_cell_dist",
                     "symphony_per_cluster_dist",
                     "predicted_doublet",]

log_message("Generando tsne plots de variables continuas...")
for color in continous_to_plot:
    fig, ax = plt.subplots()
    sc.pl.tsne(amp_random_order,
               color=color,
               ax=ax,
               show=False,
               s=size,
               alpha=1,
               frameon=True,
               )
    file_name = f"{base_file_name}_{color}"
    plt.savefig(f"{file_name}.png", bbox_inches="tight")
    plt.close()

