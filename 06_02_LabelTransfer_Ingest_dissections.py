# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 19:52:03 2024

Ingest Label Transfer

@author: Juan Andrés Tejedor Serrano
"""

# ============================================
# CARGA DE MÓDULOS
# ============================================
import argparse
import numpy as np
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
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
parser.add_argument("input_query", type=str)
parser.add_argument("input_ref", type=str)
parser.add_argument("output_query", type=str)
parser.add_argument("output_ref", type=str)
parser.add_argument("save_path", type=str)

# Extraigo los argumentos del CLI a un objeto de Python
args = parser.parse_args()

# Recupero argumentos como variables
input_query = args.input_query
input_ref = args.input_ref
output_query = args.output_query
output_ref = args.output_ref
save_path = args.save_path

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
# CONFIGURACIÓN INICIAL Y CARGA DE DATOS
# ============================================

# Establezco ajustes de Scanpy
sc.settings.verbosity = 4
sc.settings.n_jobs = 8

# Cargo datos y veo sus atributos
log_message(f"Cargando amp desde {input_query}...")
query = ad.read_h5ad(input_query)
log_message(f"Datos cargados desde {input_query}:\n")
print(query)

log_message(f"Cargando amp desde {input_ref}...")
ref = ad.read_h5ad(input_ref)
log_message(f"Datos cargados desde {input_ref}:\n")
print(ref)

"""
# ============================================
# 1. FILTRO REFERENCIA POR DISECCIONES
# ============================================

# Listo todas las disecciones que potencialmente interesan a mis datos
DMNX_dissections = [  # Medulla dissections (DMNX se halla en la medula)
    'Human MoAN', 'Human MoRF-MoEN', 'Human MoSR', 'Human IO'
]

GPI_dissections = [  # Basal forebrain dissections (GPI viene incluido explicitamente, pero también es interesante tener regiones circundantes)
    'Human CaB', 'Human GPe', 'Human GPi', 'Human NAC',
    'Human Pu', 'Human SEP', 'Human SI', 'Human Cla',
]

PMC_and_PFC_dissections = [  # Frontal + parietal cortex dissections (PMC es M1C, está incluido explicitamente, pero PFC no, por lo que incluyo áreas vecinas)
    'Human M1C', 'Human A13', 'Human A14', 'Human A25',
    'Human A32', 'Human A44-A45', 'Human A46', 'Human FI',
    'Human A5-A7', 'Human A40', 'Human A43', 'Human S1C',

]

PVC_dissections = [  # Occipital cortex dissections (PVC es V1C, está incluido explicitamente)
   'Human A19', 'Human Pro', 'Human V1C', 'Human V2'
]

# Agrupo todas
all_dissections = DMNX_dissections + GPI_dissections + PMC_and_PFC_dissections + PVC_dissections

# Filtro
ref = ref[ref.obs.roi.isin(all_dissections)]
log_message(f"Referencia filtrada por disecciones:\n{ref}...")

# ============================================
# 2. PROCESO REFERENCIA
# ============================================
# Elimino células con un abundantes o muy pocos conteos/genes
log_message("Aplicando QC...")
sc.pp.filter_cells(ref, min_counts=1500)
sc.pp.filter_cells(ref, max_counts=110000)
sc.pp.filter_cells(ref, min_genes=1100)
sc.pp.filter_cells(ref, max_genes=12500)
# Elimino células con alto contenido mitocondrial
ref = ref[ref.obs['fraction_mitochondrial']<0.02]

log_message("QC completado:\n")
print(ref)

# Almaceno conteos en objeto para luego incluir como layer adicional (la normalización afecta también a layers!)
counts = ref.X.copy()

# Aplico normalización "shifted logarithm" (size-factors + log) (library-size norm)
log_message("Normalizando...")
sc.pp.normalize_total(ref, target_sum=None, inplace=True)
log_message("Matriz de expresión post-corrección por size-factors:\n")
print(ref.X.toarray())

# Aplico ln(X+1)
sc.pp.log1p(ref)
log_message("Matriz de expresión post-normalización :\n")
print(ref.X.toarray())

# Añado layers de conteos
ref.layers["counts"] = counts

log_message("Normalización completada :\n")
print(ref)

# Selección de HVGs
log_message("Seleccionando HVGs...")
top_genes = np.floor(0.2*ref.shape[1]).astype(int)
try:
    sc.pp.highly_variable_genes(ref, layer="counts",  # El input son conteos
                                n_top_genes=top_genes,  # Top 20% HVGs
                                batch_key="donor_id",  # Considero efecto lote
                                flavor="seurat_v3",  # Creo que es más adecuado para mi caso donde hay pacientes casos y controles. Para cada paciente, este método hace un ranking de x HVGs. Tras  ello, obtiene, para cada HVG, la mediana del ranking donde ha quedado. Finalmente, obtiene un ranking de HVGs donde se ordenan por mediana (de menor (más variable) a mayor (menos variable)), seleccionándose los x primeros HVGs. Este método tiene preferencia por HVGs que son consistentes, pudiendo ignorar HVG especificos de condición. Aún así, como el nº de HVGs es bastante alto (20%) creo que no habrá problemas. 
                                span=0.3)
    log_message("Selección de HVGs completada considerando efecto lote:\n")
except Exception as e2:
    log_message(f"Error capturado al intentar seleccionar HVGs corrigiendo por efecto lote: {e2}")
    sc.pp.highly_variable_genes(ref, layer="counts",  # El input son conteos
                                n_top_genes=top_genes,  # Top 20% HVGs
                                flavor="seurat_v3",  # Creo que es más adecuado para mi caso donde hay pacientes casos y controles. Para cada paciente, este método hace un ranking de x HVGs. Tras  ello, obtiene, para cada HVG, la mediana del ranking donde ha quedado. Finalmente, obtiene un ranking de HVGs donde se ordenan por mediana (de menor (más variable) a mayor (menos variable)), seleccionándose los x primeros HVGs. Este método tiene preferencia por HVGs que son consistentes, pudiendo ignorar HVG especificos de condición. Aún así, como el nº de HVGs es bastante alto (20%) creo que no habrá problemas. 
                                span=0.3)
    log_message("Selección de HVGs completada sin considerar efecto lote:\n")
print(ref)

# PCA
log_message("Obteniendo PCA...")
sc.pp.pca(ref, n_comps=50)
log_message("PCA completado :\n")
print(ref)

# ==============================================================
# 3. CORRIJO EFECTO LOTE REFERENCIA (a nivel de participante)
# ==============================================================
# Corrigo PCs con Harmony
log_message("Ejecutando Harmony (corrigiendo PCs)...")
sc.external.pp.harmony_integrate(ref, key="donor_id")

# Reemplazo los PC de Harmony por los por defecto, aunque almaceno estos últimos
ref.obsm['X_pca_uncorrected'] = ref.obsm['X_pca']
ref.obsm['X_pca'] = ref.obsm['X_pca_harmony']

log_message("Harmony ejecutado (PCs corregidos):\n")
print(ref)

# Re-ejecuto la pipeline de clustering
# Calculo de vecinos
log_message("Obteniendo vecinos más cercanos y gráfico de vecindades...")
sc.pp.neighbors(ref, n_neighbors=30, n_pcs=50, metric="manhattan")
log_message("Cálculo de vecindades completado:\n")
print(ref)
# UMAP
log_message("Ejecutando UMAP...")
# Reemplazo UMAP corregido por el precorregido, aunque almaceno este último
sc.tl.umap(ref, min_dist=0.1)
log_message("UMAP completado:\n")
print(ref)
"""

# ============================================
# LABEL TRANSFER (vía INGEST)
# ============================================

# Antes obtengo UMAP de query
sc.tl.umap(query, min_dist=0.1)
# Almaceno el UMAP pre-ingest (el ingest lo sobreescribe por el UMAP proyectado sobre referencia)
query.obsm["X_umap_preingest"] = query.obsm["X_umap"]

# Para poder hacer la transferencia, los genes de ambos datasets deben coincidir
log_message("Seleccionando genes comunes...")
common_genes = ref.var_names.intersection(query.var_names)
query_sub = query[:, common_genes]
ref_sub = ref[:, common_genes]

log_message(f"Genes comunes seleccionados:\nQUERY:\n{query_sub}\n\nREF:\n{ref_sub}")

log_message("Transfiriendo labels de referencia a query...")
obs_to_transfr = ["ROIGroupFine",
                  "supercluster_term",
                  "cell_type",
                  "cluster_id"]
query_sub = sc.tl.ingest(query_sub, ref_sub, obs=obs_to_transfr, inplace=False)
log_message(f"Labels transferidas:\n{query_sub}")

# Transfiero las anotaciones al query original
log_message("Añadiendo anotaciones a objeto final...")

for obs in obs_to_transfr:
    query.obs[obs] = query_sub.obs[obs]
query.obs["cluster_id"] = query.obs["cluster_id"].astype("category")

log_message("Añadiendo embeddings a objeto final...")
query.obsm["X_pca_ingest"] = query_sub.obsm["X_pca"]
query.obsm["X_umap_ingest"] = query_sub.obsm["X_umap"]


sc.pp.neighbors(query, n_neighbors=30, n_pcs=50, metric="manhattan",
                use_rep="X_pca_ingest", key_added="ingest")
sc.tl.leiden(query, resolution=0.1, flavor="leidenalg",
             key_added="leiden_0_10_ingested", neighbors_key="ingest")

# ============================================
# GUARDAR EL OBJETO h5ad
# ============================================
log_message(f"Guardando el amp procesado en {output_query}...")
query.write_h5ad(output_query, compression="gzip")
log_message("Listo!")
"""
log_message(f"Guardando el ref procesado en {output_ref}...")
ref.write_h5ad(output_ref, compression="gzip")
log_message("Listo!")
"""
# ============================================
# GENERACIÓN DE PLOTS
# ============================================

# Reordeno aleatoriamente las células, evitando solapamiento de puntos.
np.random.seed(0)
random_indices = np.random.permutation(list(range(query.shape[0])))
query_random = query[random_indices, :]

# Establezco el nombre base de los archivos
plot_type = "umap_exploratory"
base_file_name = f"{save_path}amppd_{plot_type}"

categoric_to_plot = ["supercluster_term",
                     "cell_type",
                     "ROIGroupFine",
                     "leiden_1_00_corrected",
                     "leiden_0_25_corrected",
                     "leiden_0_10_corrected",
                     "leiden_0_10_ingested",
                     "set",
                     "case_control",
                     "brain_region",
                     "Phase",
                     "cluster_id"
                     ]

log_message("Generando UMAP plots de variables categóricas POST-INGEST...")
for color in categoric_to_plot:
    log_message(f"Plot: {color}")
    ncat = len(query.obs[color].cat.categories)
    fig, ax = plt.subplots()
    if ncat <= 30:
        sc.pl.embedding(query_random,
                        basis="X_umap_ingest",
                        color=color,
                        ax=ax,
                        show=False,
                        s=size,
                        alpha=1,
                        palette=sns.husl_palette(ncat),
                        frameon=True,
                        )
    else:
        sc.pl.embedding(query_random,
                        basis="X_umap_ingest",
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
    file_name = f"{base_file_name}_{color}_postingest"
    plt.savefig(f"{file_name}.png", bbox_inches="tight")
    plt.close()

log_message("Generando UMAP plots de variables categóricas PRE-INGEST...")
for color in categoric_to_plot:
    log_message(f"Plot: {color}")
    ncat = len(query.obs[color].cat.categories)
    fig, ax = plt.subplots()
    if ncat <= 30:
        sc.pl.embedding(query_random,
                        basis="X_umap_preingest",
                        color=color,
                        ax=ax,
                        show=False,
                        s=size,
                        alpha=1,
                        palette=sns.husl_palette(ncat),
                        frameon=True,
                        )
    else:
        sc.pl.embedding(query_random,
                        basis="X_umap_preingest",
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
    file_name = f"{base_file_name}_{color}_preingest"
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
                     "predicted_doublet",]

log_message("Generando UMAP plots de variables continuas POST-INGEST...")
for color in continous_to_plot:
    log_message(f"Plot: {color}")
    fig, ax = plt.subplots()
    sc.pl.embedding(query_random,
                    basis="X_umap_ingest",
                    color=color,
                    ax=ax,
                    show=False,
                    s=size,
                    alpha=1,
                    frameon=True,
                    )
    file_name = f"{base_file_name}_{color}_postingest"
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


