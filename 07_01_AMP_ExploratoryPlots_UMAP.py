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
import matplotlib
import matplotlib.pyplot as plt

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.curdir, '02-scripts')))
from my_utils import log_message
# ============================================
# PARSING DE CLI A PYTHON
# ============================================
parser = argparse.ArgumentParser(description="Plots")
parser.add_argument("input_amp", type=str)
parser.add_argument("save_path", type=str)
args = parser.parse_args()

input_amp = args.input_amp
save_path = args.save_path

# ============================================
# CONFIGURACIÓN INICIAL Y CARGA DE DATOS
# ============================================
sc.settings.verbosity = 4
sc.settings.n_jobs = 8


log_message(f"Cargando amp desde {input_amp}...")
amp = ad.read_h5ad(input_amp)
log_message(f"Datos cargados desde {input_amp}:\n")
print(amp)

# Convierto a categórica
convert = ["hoehn_yahr_stage", "path_braak_lb", "path_braak_nft"]
for var in convert:
    amp.obs[var] = amp.obs[var].astype("category")

# Reordeno aleatoriamente las células
np.random.seed(0)
random_indices = np.random.permutation(list(range(amp.shape[0])))
amp_random_order = amp[random_indices, :]


# ============================================
# CONFIGURACIÓN DE PLOTS
# ============================================
sc.settings.autoshow = False  # No mostrar gráfico en consola. Permite realizar modificaciones a los gráficos obtenidos por scanpy para después guardarlos de forma explícita.

size = 0.5  # Tamaño de puntos por defecto
highlight_size = 5  # Tamaño de puntos para los gráficos "highlight"
dpi_save = 300  # Valor por defecto de DPI
basis = "X_umap_preingest"

sc.set_figure_params(dpi=50,
                     dpi_save=dpi_save,
                     figsize=(12, 8),
                     frameon=True,
                     color_map="magma")

xkcd_colors = list(matplotlib.colors.XKCD_COLORS.values())
plt.rcParams.update({
    'axes.titlesize': 28,     # Título del gráfico
    'axes.labelsize': 20,     # Etiquetas de los ejes
    'xtick.labelsize': 16,    # Texto en el eje X
    'ytick.labelsize': 16,    # Texto en el eje Y
    'legend.fontsize': 18,    # Tamaño de la fuente de la leyenda
    'font.size': 16,           # Tamaño de fuente general
})


# ============================================
# FUNCIONES PERSONALIZADAS
# ============================================
def generate_categorical_plots(adata, save_path, basis, xkcd_colors):
    categoric_to_plot = [
        "supercluster_term",
        "cell_type",
        "ROIGroupFine",
        "leiden_1_00_corrected",
        "leiden_0_25_corrected",
        "leiden_0_10_corrected",
        "set",
        "case_control",
        "brain_region",
        "Phase",
        "hoehn_yahr_stage",
        "path_braak_lb",
        "path_braak_nft",
        "sex",
        "sample_id",
        "participant_id",
        "cohort"
    ]

    log_message("Generando umap plots de variables categóricas...")
    for color in categoric_to_plot:
        ncat = len(adata.obs[color].cat.categories)
        fig, ax = plt.subplots()
        if ncat <= 102:
            sc.pl.embedding(adata,
                            basis=basis,
                            color=color,
                            ax=ax,
                            show=False,
                            s=size,
                            alpha=1,
                            frameon=False)
        else:
            sc.pl.embedding(adata,
                            basis=basis,
                            color=color,
                            ax=ax,
                            show=False,
                            s=size,
                            alpha=1,
                            palette=xkcd_colors,
                            frameon=False)
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
        file_name = f"{save_path}amppd_umap_exploratory_{color}"
        plt.savefig(f"{file_name}.png", bbox_inches="tight")
        plt.close()

def generate_continuous_plots(adata, save_path, basis):
    continous_to_plot = [
        "log1p_n_genes_by_counts",
        "log1p_total_counts",
        "doublet_score",
        "pct_counts_in_top_20_genes",
        "percent_apop",
        "percent_dna_repair",
        "percent_ieg",
        "percent_mito",
        "percent_mito_ribo",
        "percent_oxphos",
        "predicted_doublet",
        "G2M.Score",
        "S.Score"
    ]

    log_message("Generando umap plots de variables continuas...")
    for color in continous_to_plot:
        fig, ax = plt.subplots()
        sc.pl.embedding(adata,
                        basis=basis,
                        color=color,
                        ax=ax,
                        show=False,
                        s=size,
                        alpha=1,
                        frameon=False)
        file_name = f"{save_path}amppd_umap_exploratory_{color}"
        plt.savefig(f"{file_name}.png", bbox_inches="tight")
        plt.close()

def generate_marker_gene_plots(adata, save_path, basis, marker_genes):
    log_message("Generando umap plots de genes marcadores...")
    for cell_type, markers in marker_genes.items():
        fig, axs = plt.subplots(1, len(markers), figsize=(5 * len(markers), 5))
        if len(markers) == 1:
            axs = [axs]
        for idx, marker in enumerate(markers):
            if marker in adata.var["hgnc_symbol"].values:
                sc.pl.embedding(adata,
                                basis=basis,
                                color=marker,
                                gene_symbols="hgnc_symbol",
                                ax=axs[idx],
                                show=False,
                                s=size,
                                alpha=1,
                                frameon=False)
                axs[idx].set_title(marker)
            else:
                axs[idx].text(0.5, 0.5, f"{marker} not found", ha='center', va='center')
                axs[idx].set_title(marker)
        plt.suptitle(cell_type)
        plt.tight_layout()
        file_name = f"{save_path}amppd_umap_marker_genes_{cell_type.replace(' ', '_').replace('/', '_')}"
        plt.savefig(f"{file_name}.png", bbox_inches="tight")
        plt.close()


def generate_highlight_plots(adata, save_path, conditions, basis, alpha):
    log_message("Generando umap plots...")
    for condition in conditions:
        if condition in adata.obs:
            unique_categories = adata.obs[condition].cat.categories
            log_message(f"Generando plots para la categoría {condition}")
            for category in unique_categories:
                fig, ax = plt.subplots()
                sc.pl.umap(adata,
                           color=[condition],
                           groups=[category],
                           ax=ax,
                           show=False,
                           s=highlight_size,
                           alpha=alpha)
                plot_type = f"umap_highlight_{condition}_{category}"
                file_name = f"{save_path}amppd_{plot_type}"
                plt.savefig(f"{file_name}.png", bbox_inches="tight")
                plt.close()


generate_categorical_plots(amp_random_order, save_path, basis, xkcd_colors)


generate_continuous_plots(amp_random_order, save_path, basis)


marker_genes = {
    "Glutamatergic Neurons (GLU)": ["SLC17A7", "GRIN1", "CAMK2A", "VGLUT1"],
    "GABAergic Neurons (GABA)": ["GAD1", "GAD2", "SLC32A1", "VGAT"],
    "Interneurons (IN)": ["SST", "PVALB", "VIP"],
    "Astrocytes (AST)": ["GFAP", "AQP4", "SLC1A3", "ALDH1L1"],
    "Oligodendrocytes (OL)": ["MBP", "MOG", "OLIG1", "PLP1"],
    "Microglia (MG)": ["AIF1", "CX3CR1", "TMEM119", "P2RY12"],
    "Ependymal Cells (EP)": ["S100B", "VIM", "FOXJ1", "TTR"],
    "Vascular Cells (VC)": ["CLDN5", "PECAM1", "VWF", "CD34"],
    "Perivascular Cells (PV)": ["PDGFRB", "RGS5", "ANPEP", "ACTA2"],
    "Dopaminergic Neurons (DA)": ["TH", "DAT", "Nurr1", "SLC18A2", "SLC6A3", "NR4A2"],
    "Serotonergic Neurons (5HT)": ["SLC6A4", "TPH2", "HTR2A"],
    "Cholinergic Neurons (ACh)": ["CHAT", "VAChT", "ChAT", "SLC18A3"],
    "Noradrenergic Neurons (NA)": ["DBH", "NET", "TH", "SLC6A2"],
    "Histaminergic Neurons (HA)": ["HDC", "HRH3", "HNMT"],
    "Schwann Cells (SC)": ["MPZ", "PMP22", "SOX10", "S100B"],
    "Deep Layer CT/6b Neurons (DLCT6B)": ["TLE4", "FOXP2", "BCL11B"],
    "Deep Layer Near Projecting Neurons (DLNP)": ["FEZF2", "CRYM", "TBR1"],
    "Hippocampal CA4 Neurons (CA4)": ["PROX1", "C1QL2"],
    "Hippocampal CA1-3 Neurons (CA1-3)": ["PCP4", "RGS14"],
    "Amygdala Excitatory Neurons (AMY)": ["SATB2", "TBR1"],
    "Deep Layer Intratelencephalic Neurons (DLIT)": ["CUX2", "NR4A2"],
    "MGE Interneurons (MGE-IN)": ["LHX6", "NKX2-1", "SST"],
    "LAMP5-LHX6 & Chandelier Neurons (CHN)": ["LAMP5", "LHX6", "PVALB"],
    "Upper Rhombic Lip Cells (URLC)": ["ATOH1", "BARHL1"],
    "Upper Layer Intratelencephalic Neurons (ULIT)": ["CUX1", "RORB", "SATB2"],
    "Thalamic Excitatory Neurons (THAL)": ["SLC17A6", "SLC17A7"],
    "CGE Interneurons (CGE-IN)": ["VIP", "CRH", "SP8"],
    "Midbrain Inhibitory Neurons (MB-IN)": ["GAD1", "GAD2", "TH"],
    "Medium Spiny Neurons (MSN)": ["DARPP-32", "PPP1R1B", "DRD1", "DRD2"],
    "Eccentric Medium Spiny Neurons (E-MSN)": ["FOXP1", "EBF1"],
    "Oligodendrocyte Precursor Cells (OPC)": ["PDGFRA", "NG2", "CSPG4", "OLIG2"],
    "Fibroblasts (FIB)": ["COL1A1", "PDGFRA", "FSP1"],
    "Choroid Plexus (CP)": ["TTR", "KCNJ13", "AQP1", "CST3"],
    "Parkinson Disease (PD)": ["SNCA", "PARK7", "LRRK2", "GBA", "TH"]
}
generate_marker_gene_plots(amp_random_order, save_path, basis, marker_genes)


conditions = ["supercluster_term",
              "cell_type",
              "leiden_0_10_corrected",
              "set",
              "participant_id",
              "brain_region",
              ]
generate_highlight_plots(amp_random_order, save_path, conditions, basis, alpha=1)
