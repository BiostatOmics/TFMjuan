# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 16:02:25 2024

Batch Correction and Visualization

@author: Juan Andrés Tejedor Serrano
"""

# ============================================
# CARGA DE MÓDULOS
# ============================================
import argparse
import numpy as np
import anndata as ad
import scanpy as sc
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
parser = argparse.ArgumentParser(description="Processing")

# Establezco argumentos a parsear
parser.add_argument("input_amp", type=str)
parser.add_argument("output_amp", type=str)
parser.add_argument("save_path", type=str)

# Extraigo los argumentos del CLI a un objeto de Python
args = parser.parse_args()

# Recupero argumentos como variables
input_amp = args.input_amp
output_amp = args.output_amp
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
sc.set_figure_params(dpi=80,
                     dpi_save=300,  # DPI de 300 para alta calidad
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

# Reordeno aleatoreamiente las células, evitando solapamiento de puntos.
np.random.seed(0)
random_indices = np.random.permutation(list(range(amp.shape[0])))
amp_random_order = amp[random_indices, :]


# ============================================
# VISUALIZACIÓN PRE-CORRECCIÓN
# ============================================


log_message("Generando tSNE plots...")

# Establezco el nombre base de los archivos
plot_type = "tsne_precorrected"
base_file_name = f"{save_path}amppd_{plot_type}"

to_plot = ["leiden_0_10", "participant_id", "sample_id",
           "set", "case_control", "brain_region"]

for color in to_plot:
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



# ============================================
# CORRECCIÓN CON HARMONY
# ============================================

# Corrigo PCs con Harmony
log_message("Ejecutando Harmony (corrigiendo PCs)...")
sc.external.pp.harmony_integrate(amp, key="participant_id")

# Reemplazo los PC de Harmony por los por defecto, aunque almaceno estos últimos
amp.obsm['X_pca_uncorrected'] = amp.obsm['X_pca']
amp.obsm['X_pca'] = amp.obsm['X_pca_harmony']

log_message("Harmony ejecutado (PCs corregidos):\n")
print(amp)

# Re-ejecuto la pipeline de clustering
# Calculo de vecinos
log_message("Obteniendo vecinos más cercanos y gráfico de vecindades...")
sc.pp.neighbors(amp, n_neighbors=30, n_pcs=50, metric="manhattan")
log_message("Cálculo de vecindades completado:\n")
print(amp)
# Calculo de clústeres
log_message("Ejecutando Leiden clustering a res = 0.1...")
sc.tl.leiden(amp, resolution=0.1, flavor="leidenalg", key_added="leiden_0_10_corrected")
log_message("Ejecutando Leiden clustering a res = 0.25...")
sc.tl.leiden(amp, resolution=0.25, flavor="leidenalg", key_added="leiden_0_25_corrected")
log_message("Ejecutando Leiden clustering a res = 1...")
sc.tl.leiden(amp, resolution=1, flavor="leidenalg", key_added="leiden_1_00_corrected")
log_message("Clustering completado:\n")
print(amp)
# tSNE
log_message("Ejecutando tSNE...")
amp.obsm["X_tsne_uncorrected"] = amp.obsm["X_tsne"]
sc.tl.tsne(amp, n_pcs=50, perplexity=500, early_exaggeration=20)
log_message("tSNE completado:\n")
print(amp)

# ============================================
# GUARDAR EL OBJETO h5ad
# ============================================
log_message(f"Guardando el amp procesado en {output_amp}...")
amp.write_h5ad(output_amp, compression="gzip")
log_message("Listo!")


# ============================================
# VISUALIZACIÓN POST-CORRECCIÓN
# ============================================

# Reordeno aleatoreamiente las células, evitando solapamiento de puntos.
np.random.seed(0)
random_indices = np.random.permutation(list(range(amp.shape[0])))
amp_random_order = amp[random_indices, :]

# Establezco el nombre base de los archivos
plot_type = "tsne_postcorrected"
base_file_name = f"{save_path}amppd_{plot_type}"

to_plot = ["leiden_0_10", "participant_id", "sample_id",
           "set", "case_control", "brain_region"]

for color in to_plot:
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


