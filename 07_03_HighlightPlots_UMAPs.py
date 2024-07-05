# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 12:19:28 2024

Highlight Draft Plots for Exploration

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
dpi_save = 500
sc.set_figure_params(dpi=80,
                     dpi_save=dpi_save,  # DPI de 300 para alta calidad
                     figsize=(12, 8),  # Tamaño de 14 pulgadas x 10 pulgadas
                     frameon=True,
                     color_map="magma")  # Mostrar marco (ejes)

xkcd_colors = list(matplotlib.colors.XKCD_COLORS.values())
plt.rcParams.update({
    'axes.titlesize': 28,     # Título del gráfico
    'axes.labelsize': 20,     # Etiquetas de los ejes
    'xtick.labelsize': 16,    # Texto en el eje X
    'ytick.labelsize': 16,    # Texto en el eje Y
    'legend.fontsize': 18,    # Tamaño de la fuente de la leyenda
    'font.size': 16,           # Tamaño de fuente general
})
size = 0.5


beauty_colors_ext = [
    "#EF476F",
    "#118AB2",
    "#06D6A0",
    "#FFD166",
    "#073B4C",
    "#FF5DF3",
    "#ff8e5d",
    "#8eba00",
    "#79d7ff",
    "#D33322",
    "#964d51",
    "#C46904",
    "#5d8e75",
    "#565656",
    "#790000",
    "#A40052",
    "#008200",
    "#975EFE"
]


not_beauty_colors = glasbey.extend_palette(palette=["#2274A5", "#F75C03", "#F1C40F", "#D90368", "#00CC66"], palette_size=40)
# Crear una paleta de colores personalizada
cell_types = amp.obs['annotations'].cat.categories
num_colors = len(cell_types)

if len(beauty_colors_ext) < num_colors:
    beauty_colors_ext = glasbey.extend_palette(palette=beauty_colors_ext, palette_size=num_colors)

# Asignar colores a las categorías (a las celulas mas abundantes se le asignan los primeros colores)
cell_types_by_abundance = amp.obs['annotations'].cat.reorder_categories(list(amp.obs['annotations'].value_counts().index)).cat.categories
color_map = {cell_type: beauty_colors_ext[i] for i, cell_type in enumerate(cell_types_by_abundance)}

# ============================================
# GENERACIÓN DE PLOTS
# ============================================

# Reordeno aleatoriamente las células, evitando solapamiento de puntos.
np.random.seed(0)
random_indices = np.random.permutation(list(range(amp.shape[0])))
amp_random_order = amp[random_indices, :]

def plot_and_save_umap(data, save_path, condition, category,
                       plot_type, dpi, s, alpha):
    base_file_name = f"{save_path}amppd_{plot_type}"
    fig, ax = plt.subplots()
    sc.pl.embedding(data,
                    basis = "X_umap_preingest",
               color=[condition],
               groups=[category],
               ax=ax,
               show=False,
               s=s,
               alpha=alpha, 
               frameon=False)
    file_name = f"{base_file_name}_{category}"
    plt.savefig(f"{file_name}.png", dpi=dpi, bbox_inches="tight")
    plt.close()
    log_message(f"Plot generado para {condition} - {category}")


def generate_highlight_umap_plots(amp, save_path, conditions,
                                  dpi=500, s=0.5, alpha=1):
    log_message("Generando umap plots...")

    for condition in conditions:
        if condition in amp.obs:
            unique_categories = amp.obs[condition].cat.categories
            plot_type = f"umap_highlight_{condition}"
            log_message(f"Generando plots para la categoría {condition}")
            for category in unique_categories:
                plot_and_save_umap(amp,save_path, condition,
                                   category, plot_type, dpi, s, alpha)


conditions = ["supercluster_term",
              #"leiden_0_10_corrected",
              #"set",
              #"participant_id",
              "brain_region",
              "case_control"
              ]

generate_highlight_umap_plots(amp, save_path=save_path,
                              conditions=conditions, dpi=500, s=2.5, alpha=1)
