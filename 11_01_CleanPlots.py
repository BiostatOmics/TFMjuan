# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 22:28:18 2024

Clean UMAPs and more!

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
import glasbey
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

# Filtro y elimino las células Unknown
amp = amp[amp.obs["annotations"] != "Unknown"].copy()

# Convierto a categórica
convert = ["hoehn_yahr_stage", "path_braak_lb", "path_braak_nft"]
for var in convert:
    amp.obs[var] = amp.obs[var].astype("category")


# Obtengo UMAP pre corrección efecto lote
sc.pp.neighbors(amp, n_neighbors=30, n_pcs=50,
                      use_rep="X_pca_uncorrected", metric="manhattan",
                      key_added="preharmony")
temp = sc.tl.umap(amp, min_dist=0.1, copy=True, neighbors_key="preharmony")
amp.obsm["X_umap_preharmony"] = temp.obsm["X_umap"]

basis_pre = "X_umap_preharmony"
basis_post = "X_umap_preingest"


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


# ============================================
# GENERACIÓN DE PLOTS
# ============================================


# ============================================
# SETUP
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


# Reordeno aleatoriamente las células
np.random.seed(0)
random_indices = np.random.permutation(list(range(amp.shape[0])))
amp = amp[random_indices, :]

# ============================================
# PLOTS
# ============================================
# CELL TYPES
basis = "X_umap_preingest"
color = "annotations"
fig, ax = plt.subplots()
sc.pl.embedding(amp,
                basis=basis,
                color=color,
                ax=ax,
                show=False,
                s=size,
                alpha=1,
                frameon=False,
                palette=color_map,
                title="")
file_name = f"{save_path}amppd_umap_{color}_{basis}"
plt.savefig(f"{file_name}.png", bbox_inches="tight")

# ============================================
# LEIDEN 0.1
basis = "X_umap_preingest"
color = "leiden_0_10_corrected"
fig, ax = plt.subplots()
sc.pl.embedding(amp,
                basis=basis,
                color=color,
                ax=ax,
                show=False,
                legend_loc="on data",
                legend_fontoutline=2,
                s=size,
                alpha=1,
                frameon=False,
                palette=not_beauty_colors,
                title="")
file_name = f"{save_path}amppd_umap_{color}_{basis}"
plt.savefig(f"{file_name}.png", bbox_inches="tight")
# LEIDEN 0.25
basis = "X_umap_preingest"
color = "leiden_0_25_corrected"
fig, ax = plt.subplots()
sc.pl.embedding(amp,
                basis=basis,
                color=color,
                ax=ax,
                show=False,
                legend_loc="on data",
                legend_fontoutline=2,
                s=size,
                alpha=1,
                frameon=False,
                palette=not_beauty_colors,
                title="")
file_name = f"{save_path}amppd_umap_{color}_{basis}"
plt.savefig(f"{file_name}.png", bbox_inches="tight")
# LEIDEN 1.0
basis = "X_umap_preingest"
color = "leiden_1_00_corrected"
fig, ax = plt.subplots()
sc.pl.embedding(amp,
                basis=basis,
                color=color,
                ax=ax,
                show=False,
                legend_loc="on data",
                legend_fontoutline=2,
                s=size,
                alpha=1,
                frameon=False,
                palette=not_beauty_colors,
                title="")
file_name = f"{save_path}amppd_umap_{color}_{basis}"
plt.savefig(f"{file_name}.png", bbox_inches="tight")


# ============================================
# BRAIN REGION
basis = "X_umap_preingest"
color = "brain_region"
fig, ax = plt.subplots()
sc.pl.embedding(amp,
                basis=basis,
                color=color,
                ax=ax,
                show=False,
                s=size,
                alpha=1,
                frameon=False,
                palette=beauty_colors_ext,
                title="")
file_name = f"{save_path}amppd_umap_{color}_{basis}"
plt.savefig(f"{file_name}.png", bbox_inches="tight")
# ============================================
# CASE CONTROL
basis = "X_umap_preingest"
color = "case_control"
fig, ax = plt.subplots()
sc.pl.embedding(amp,
                basis=basis,
                color=color,
                ax=ax,
                show=False,
                s=size,
                alpha=1,
                frameon=False,
                palette=beauty_colors_ext,
                title="")
file_name = f"{save_path}amppd_umap_{color}_{basis}"
plt.savefig(f"{file_name}.png", bbox_inches="tight")


# ============================================
# SEX
basis = "X_umap_preingest"
color = "sex"
fig, ax = plt.subplots()
sc.pl.embedding(amp,
                basis=basis,
                color=color,
                ax=ax,
                show=False,
                s=size,
                alpha=1,
                frameon=False,
                palette=beauty_colors_ext,
                title="")
file_name = f"{save_path}amppd_umap_{color}_{basis}"
plt.savefig(f"{file_name}.png", bbox_inches="tight")


# ============================================
# DOUBLETS
basis = "X_umap_preingest"
color = "predicted_doublet"
fig, ax = plt.subplots()
sc.pl.embedding(amp,
                basis=basis,
                color=color,
                ax=ax,
                show=False,
                s=size,
                alpha=1,
                frameon=False,
                palette=beauty_colors_ext,
                title="")
file_name = f"{save_path}amppd_umap_{color}_{basis}"
plt.savefig(f"{file_name}.png", bbox_inches="tight")


# ============================================
# BRAAK
basis = "X_umap_preingest"
color = "path_braak_lb"
fig, ax = plt.subplots()
sc.pl.embedding(amp,
                basis=basis,
                color=color,
                ax=ax,
                show=False,
                s=size,
                alpha=1,
                frameon=False,
                palette=beauty_colors_ext,
                title="")
file_name = f"{save_path}amppd_umap_{color}_{basis}"
plt.savefig(f"{file_name}.png", bbox_inches="tight")

# ============================================
# BATCH EFFECT: PARTICIPANT, SAMPLE, SET
basis = "X_umap_preingest"
color = "sample_id"
fig, ax = plt.subplots()
sc.pl.embedding(amp,
                basis=basis,
                color=color,
                ax=ax,
                show=False,
                s=size,
                alpha=1,
                frameon=False,
                palette=xkcd_colors,
                title="")
plt.legend(ncol=6, bbox_to_anchor=(1, 1), fontsize=6, frameon=False)
file_name = f"{save_path}amppd_umap_{color}_{basis}"
plt.savefig(f"{file_name}.png", bbox_inches="tight")


color = "participant_id"
fig, ax = plt.subplots()
sc.pl.embedding(amp,
                basis=basis,
                color=color,
                ax=ax,
                show=False,
                s=size,
                alpha=1,
                frameon=False,
                palette=xkcd_colors,
                title="")
plt.legend(ncol=4, bbox_to_anchor=(1, 1), frameon=False)
file_name = f"{save_path}amppd_umap_{color}_{basis}"
plt.savefig(f"{file_name}.png", bbox_inches="tight")


color = "set"
fig, ax = plt.subplots()
sc.pl.embedding(amp,
                basis=basis,
                color=color,
                ax=ax,
                show=False,
                s=size,
                alpha=1,
                frameon=False,
                palette=xkcd_colors,
                title="")
plt.legend(ncol=3, bbox_to_anchor=(1, 1), frameon=False)
file_name = f"{save_path}amppd_umap_{color}_{basis}"
plt.savefig(f"{file_name}.png", bbox_inches="tight")


color = "cohort"
fig, ax = plt.subplots()
sc.pl.embedding(amp,
                basis=basis,
                color=color,
                ax=ax,
                show=False,
                s=size,
                alpha=1,
                frameon=False,
                palette=beauty_colors_ext,
                title="")
file_name = f"{save_path}amppd_umap_{color}_{basis}"
plt.savefig(f"{file_name}.png", bbox_inches="tight")


basis = "X_umap_preharmony"
color = "sample_id"
fig, ax = plt.subplots()
sc.pl.embedding(amp,
                basis=basis,
                color=color,
                ax=ax,
                show=False,
                s=size,
                alpha=1,
                frameon=False,
                palette=xkcd_colors,
                title="")
plt.legend(ncol=6, bbox_to_anchor=(1, 1), fontsize=6, frameon=False)
file_name = f"{save_path}amppd_umap_{color}_{basis}"
plt.savefig(f"{file_name}.png", bbox_inches="tight")


color = "participant_id"
fig, ax = plt.subplots()
sc.pl.embedding(amp,
                basis=basis,
                color=color,
                ax=ax,
                show=False,
                s=size,
                alpha=1,
                frameon=False,
                palette=xkcd_colors,
                title="")
plt.legend(ncol=4, bbox_to_anchor=(1, 1), frameon=False)
file_name = f"{save_path}amppd_umap_{color}_{basis}"
plt.savefig(f"{file_name}.png", bbox_inches="tight")


color = "set"
fig, ax = plt.subplots()
sc.pl.embedding(amp,
                basis=basis,
                color=color,
                ax=ax,
                show=False,
                s=size,
                alpha=1,
                frameon=False,
                palette=xkcd_colors,
                title="")
plt.legend(ncol=3, bbox_to_anchor=(1, 1), frameon=False)
file_name = f"{save_path}amppd_umap_{color}_{basis}"
plt.savefig(f"{file_name}.png", bbox_inches="tight")


color = "cohort"
fig, ax = plt.subplots()
sc.pl.embedding(amp,
                basis=basis,
                color=color,
                ax=ax,
                show=False,
                s=size,
                alpha=1,
                frameon=False,
                palette=beauty_colors_ext,
                title="")
file_name = f"{save_path}amppd_umap_{color}_{basis}"
plt.savefig(f"{file_name}.png", bbox_inches="tight")


# ============================================
# COUNTS & GENES
basis = "X_umap_preingest"
color = "total_counts"
fig, ax = plt.subplots()
sc.pl.embedding(amp,
                basis=basis,
                color=color,
                ax=ax,
                show=False,
                s=size,
                alpha=1,
                frameon=False,
                title="")
plt.legend(ncol=6, bbox_to_anchor=(1, 1), fontsize=6, frameon=False)
file_name = f"{save_path}amppd_umap_{color}_{basis}"
plt.savefig(f"{file_name}.png", bbox_inches="tight")

basis = "X_umap_preingest"
color = "n_genes_by_counts"
fig, ax = plt.subplots()
sc.pl.embedding(amp,
                basis=basis,
                color=color,
                ax=ax,
                show=False,
                s=size,
                alpha=1,
                frameon=False,
                title="")
plt.legend(ncol=6, bbox_to_anchor=(1, 1), fontsize=6, frameon=False)
file_name = f"{save_path}amppd_umap_{color}_{basis}"
plt.savefig(f"{file_name}.png", bbox_inches="tight")


basis = "X_umap_preingest"
color = "log1p_total_counts"
fig, ax = plt.subplots()
sc.pl.embedding(amp,
                basis=basis,
                color=color,
                ax=ax,
                show=False,
                s=size,
                alpha=1,
                frameon=False,
                title="")
plt.legend(ncol=6, bbox_to_anchor=(1, 1), fontsize=6, frameon=False)
file_name = f"{save_path}amppd_umap_{color}_{basis}"
plt.savefig(f"{file_name}.png", bbox_inches="tight")

basis = "X_umap_preingest"
color = "log1p_n_genes_by_counts"
fig, ax = plt.subplots()
sc.pl.embedding(amp,
                basis=basis,
                color=color,
                ax=ax,
                show=False,
                s=size,
                alpha=1,
                frameon=False,
                title="")
plt.legend(ncol=6, bbox_to_anchor=(1, 1), fontsize=6, frameon=False)
file_name = f"{save_path}amppd_umap_{color}_{basis}"
plt.savefig(f"{file_name}.png", bbox_inches="tight")


# ============================================
# MARKER GENES
marker_genes = [
    "SLC17A7", "GRIN1", "CAMK2A", "VGLUT1",
    "GAD1", "GAD2", "SLC32A1", "VGAT",
    "SST", "PVALB", "VIP",
    "GFAP", "AQP4", "SLC1A3", "ALDH1L1",
    "MBP", "MOG", "OLIG1", "PLP1",
    "AIF1", "CX3CR1", "TMEM119", "P2RY12",
    "S100B", "VIM", "FOXJ1", "TTR",
    "CLDN5", "PECAM1", "VWF", "CD34",
    "PDGFRB", "RGS5", "ANPEP", "ACTA2",
    "TH", "DAT", "Nurr1", "SLC18A2", "SLC6A3", "NR4A2",
    "SLC6A4", "TPH2", "HTR2A",
    "CHAT", "VAChT", "ChAT", "SLC18A3",
    "DBH", "NET", "SLC6A2",
    "HDC", "HRH3", "HNMT",
    "MPZ", "PMP22", "SOX10",
    "TLE4", "FOXP2", "BCL11B",
    "FEZF2", "CRYM", "TBR1",
    "PROX1", "C1QL2",
    "PCP4", "RGS14",
    "SATB2", "CUX2",
    "LHX6", "NKX2-1",
    "LAMP5", "ATOH1", "BARHL1",
    "CUX1", "RORB",
    "SLC17A6", "CRH", "SP8",
    "DARPP-32", "PPP1R1B", "DRD1", "DRD2",
    "FOXP1", "EBF1",
    "PDGFRA", "NG2", "CSPG4", "OLIG2",
    "COL1A1", "FSP1",
    "KCNJ13", "AQP1", "CST3",
    "SNCA", "PARK7", "LRRK2", "GBA",
    # Additional marker genes
    "MAP2", "NEUN", "GLAST", "IBA1", "CNP", "MAG", "CD31", "VE-Cadherin", "Nestin", "SOX2"
]

basis = "X_umap_preingest"
for marker in marker_genes:
    fig, ax = plt.subplots()
    if marker in amp.var["hgnc_symbol"].values:
        sc.pl.embedding(amp,
                        basis=basis,
                        color=marker,
                        gene_symbols="hgnc_symbol",
                        ax=ax,
                        show=False,
                        s=size,
                        alpha=1,
                        frameon=False)
        ax.set_title(marker)
    else:
        ax.text(0.5, 0.5, f"{marker} not found", ha='center', va='center')
        ax.set_title(marker)
    plt.tight_layout()
    file_name = f"{save_path}umap_marker_genes_{marker.replace(' ', '_').replace('/', '_')}"
    plt.savefig(f"{file_name}.png", bbox_inches="tight")
    plt.close()
