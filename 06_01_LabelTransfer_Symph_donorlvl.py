# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 23:42:30 2024 (V1)
Created on Fri Jun 28 12:55:14 2024 (V2)

Symphony Automatic Annotation V2

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
parser = argparse.ArgumentParser(description="Processing")

# Establezco argumentos a parsear
parser.add_argument("input_query", type=str)
parser.add_argument("input_ref", type=str)
parser.add_argument("output_query", type=str)
parser.add_argument("save_path", type=str)
# Extraigo los argumentos del CLI a un objeto de Python
args = parser.parse_args()

# Recupero argumentos como variables
input_query = args.input_query
input_ref = args.input_ref
output_query = args.output_query
save_path = args.save_path

# ============================================
# CONFIGURACIÓN INICIAL Y CARGA DE DATOS
# ============================================

# Establezco ajustes de Scanpy
sc.settings.verbosity = 4
sc.settings.n_jobs = 8


# Cargo datos y veo sus atributos
log_message(f"Cargando query desde {input_query}...")
query = ad.read_h5ad(input_query)
log_message(f"Datos cargados desde {input_query}:\n")
print(query)

log_message(f"Cargando ref desde {input_ref}...")
ref = ad.read_h5ad(input_ref)
log_message(f"Datos cargados desde {input_ref}:\n")
print(ref)

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


# ============================================
# ANOTACIÓN AUTOMÁTICA
# ============================================


log_message("### 1. Comenzando anotación automática... ###")

log_message("2. Preparando referencia...")

log_message("2.1 Eliminando efecto lote con Harmony...")
sp.pp.harmony_integrate(ref, key="donor_id")

log_message("2.2 Recalculando vecinos...")
sc.pp.neighbors(ref, use_rep="X_pca_harmony", n_neighbors=30,
                n_pcs=50, metric="manhattan")

log_message("2.3 Obteniendo UMAP...")
sc.tl.umap(ref, min_dist=0.35)
log_message(f"Referencia lista:\n{ref}")

log_message("3. Mapeando embedding (PCA corregido) de query a referencia...")
# Añado columna para que su nombre coincida con el esperado por symphony
ref.var["mean"] = ref.var["means"]
ref.var["std"] = np.sqrt(ref.var["variances"])
sp.tl.map_embedding(query, ref, key="participant_id")
log_message(f"Embedding mapeado:\nQUERY:\n{query}\n\nREF:\n{ref}")


log_message("4. Mapeando embedding (UMAP) de query a referencia...")
sp.tl.ingest(query, ref)

log_message("5. Transfiriendo labels de referencia a query...")
sp.tl.transfer_labels_kNN(query, ref,
                          ref_labels=["supercluster_term",
                                      "cell_type",
                                      "ROIGroupFine"])
log_message(f"Labels transferidas:\n{query}")

log_message("6. Calculando scores de query a nivel de célula...")
sp.tl.per_cell_confidence(query, ref)
log_message(f"Scores calculados:\n{query}")


log_message("7. Calculando scores de query a nivel de clúster...")

log_message("7.1 Obteniendo vecinos de query integrado (a nivel de harmony)...")
sc.pp.neighbors(query, use_rep="X_pca_harmony",
                n_neighbors=30, n_pcs=50, metric="manhattan", key_added="harm")

log_message("7.2 Obteniendo clústers de query integrado (a nivel de harmony)...")
sc.tl.leiden(query, resolution=0.25, flavor="leidenalg",
             key_added="leiden_0_25_harm", neighbors_key="harm")


log_message("7.3 Obtención de scores (a nivel de harmony)...")
sp.tl.per_cluster_confidence(query, ref,
                             cluster_key="leiden_0_25_harm",
                             )


log_message(f"Scores calculados:\n{query}")


# ============================================
# GUARDAR EL OBJETO h5ad
# ============================================
log_message(f"Guardando el query procesado en {output_query}...")
# Elimino los clusters scores de uns (da error), pues ya están en obs.
del query.uns["symphony_per_cluster_dist"]
query.write_h5ad(output_query, compression="gzip")
log_message("Listo!")


# ============================================
# GENERACIÓN DE PLOTS
# ============================================

# Reordeno aleatoriamente las células, evitando solapamiento de puntos.
np.random.seed(0)
random_indices = np.random.permutation(list(range(query.shape[0])))
query_random = query[random_indices, :]

# Establezco el nombre base de los archivos
plot_type = "umap_exploratory_symphony"
base_file_name = f"{save_path}amppd_{plot_type}"

categoric_to_plot = ["supercluster_term",
                     "cell_type",
                     "ROIGroupFine",
                     "leiden_0_25",  # este es el original, sin corregir ni por lote ni integrado
                     "leiden_0_25_harm",  # este tal vez sea el que quiera aunque creo es el siguiente
                     "set",
                     "case_control",
                     "brain_region",
                     "Phase",
                     "sample_id",
                     "participant_id"
                     ]

log_message("Generando UMAP plots de variables categóricas...")
for color in categoric_to_plot:
    log_message(f"Plot: {color}")
    ncat = len(query.obs[color].cat.categories)
    fig, ax = plt.subplots()
    if ncat <= 30:
        sc.pl.umap(query_random,
                   color=color,
                   ax=ax,
                   show=False,
                   s=size,
                   alpha=1,
                   palette=sns.husl_palette(ncat),
                   frameon=True,
                   )
    else:
        sc.pl.umap(query_random,
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
                     "predicted_doublet",
                     "symphony_per_cell_dist",
                     "symphony_per_cluster_dist"
                     ]

log_message("Generando UMAP plots de variables continuas...")
for color in continous_to_plot:
    log_message(f"Plot: {color}")
    fig, ax = plt.subplots()
    sc.pl.umap(query_random,
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


# Reordeno aleatoriamente las células, evitando solapamiento de puntos.
np.random.seed(0)
random_indices = np.random.permutation(list(range(ref.shape[0])))
ref_random = ref[random_indices, :]
log_message("Generando UMAP plot de referencia...")
base_file_name = f"{save_path}ref_{plot_type}"
fig, ax = plt.subplots()
sc.pl.umap(ref_random,
           color="supercluster_term",
           ax=ax,
           show=False,
           s=size,
           alpha=1,
           frameon=True,
           )
file_name = f"{base_file_name}_supercluster_term"
plt.savefig(f"{file_name}.png", bbox_inches="tight")
plt.close()
