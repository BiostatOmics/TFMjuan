# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 00:17:44 2024

Gene Metadata Transfer, Normalization, HVG Selection, PCA, Clustering and tSNE

@author: Juan Andrés Tejedor Serrano
"""

# ============================================
# CARGA DE MÓDULOS
# ============================================
import argparse
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt

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
parser.add_argument("input_gene_metadata", type=str)
parser.add_argument("output_amp", type=str)

# Extraigo los argumentos del CLI a un objeto de Python
args = parser.parse_args()

# Recupero argumentos como variables
input_amp = args.input_amp
input_gene_metadata = args.input_gene_metadata
output_amp = args.output_amp

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

log_message(f"Cargando metadatos génicos desde {input_gene_metadata}...")
gene_metadata = pd.read_csv(input_gene_metadata, sep="\t") # 04-meta_data/02-gene_ids.tsv
log_message(f"Datos cargados desde {input_gene_metadata}:\n")
print(gene_metadata)


# ============================================
# TRANSFERENCIA DE METADATOS GÉNICOS
# ============================================

log_message("Añadiendo metadatos génicos...")
# Establezco el ENSEMBL ID como nombre de filas (index)
gene_metadata.index = gene_metadata['original_ensembl_id']
# Extraigo únicamente el GENE NAME, su DESCRIPCIÓN y el ENSEMBL ID (index)
gene_metadata = gene_metadata[['hgnc_symbol', 'description']]
# Transfiero GENE NAME y DESCRIPCIÓN a metadatos de genes, usando como anclaje ENSEMBL ID
amp.var = amp.var.merge(gene_metadata, "inner", validate="one_to_one",
                        left_index=True, right_index=True)
log_message("Metadatos génicos transferidos:\n")
print(amp.var)


# ============================================
# DETECCIÓN DE DOUBLETS
# ============================================

log_message("Detectando doublets usando scrublet...")
# Intentar corregir efecto lote por participante, de lo contrario, no corregir
try:
    amp = sc.pp.scrublet(amp, n_prin_comps=50, batch_key="participant_id",
                         copy=True, verbose=True)
    log_message("Detección completada corregiendo efecto lote por participante")
except Exception as e1:
    log_message(f"Error capturado al intentar corregir efecto lote en scrublet por participante: {e1}")
    amp = sc.pp.scrublet(amp, n_prin_comps=50,
                         copy=True, verbose=True)
    log_message("Detección completada sin corregir efecto lote")

log_message("Detección de doublets completada:\n")
print(amp)

# ============================================
# CÁLCULO DE MÉTRICAS DE QC
# ============================================

# Cálculo de otras métricas como top_genes
log_message("Calculando métricas de QC...")
sc.pp.calculate_qc_metrics(amp, percent_top=[20], inplace=True)
log_message("Métricas calculadas:\n")
print(amp)

# ============================================
# NORMALIZACIÓN
# ============================================
log_message("Matriz de expresión pre-normalización:\n")
print(amp.X.toarray())

# Almaceno conteos en objeto para luego incluir como layer adicional (la normalización afecta también a layers!)
counts = amp.X.copy()

# Aplico normalización "shifted logarithm" (size-factors + log) (library-size norm)
log_message("Normalizando...")
sc.pp.normalize_total(amp, target_sum=None, inplace=True)
log_message("Matriz de expresión post-corrección por size-factors:\n")
print(amp.X.toarray())

# Almaceno datos corregidos por el size-factor (para cada célula, dividir los conteos de cada gen por el número de conteos total de dicha célula)
cpm = amp.X.copy()

# Aplico ln(X+1)
sc.pp.log1p(amp)
log_message("Matriz de expresión post-normalización :\n")
print(amp.X.toarray())

# Añado layers de conteos y cpm
amp.layers["counts"] = counts
amp.layers["cpm"] = cpm

log_message("Normalización completada :\n")
print(amp)


# ============================================
# SELECCIÓN DE HVGs
# ============================================

# Selección de HVGs
log_message("Seleccionando HVGs...")
top_genes = np.floor(0.2*amp.shape[1]).astype(int)
try:
    sc.pp.highly_variable_genes(amp, layer="counts",  # El input son conteos
                                n_top_genes=top_genes,  # Top 20% HVGs
                                batch_key="participant_id",  # Considero efecto lote
                                flavor="seurat_v3",  # Creo que es más adecuado para mi caso donde hay pacientes casos y controles. Para cada paciente, este método hace un ranking de x HVGs. Tras  ello, obtiene, para cada HVG, la mediana del ranking donde ha quedado. Finalmente, obtiene un ranking de HVGs donde se ordenan por mediana (de menor (más variable) a mayor (menos variable)), seleccionándose los x primeros HVGs. Este método tiene preferencia por HVGs que son consistentes, pudiendo ignorar HVG especificos de condición. Aún así, como el nº de HVGs es bastante alto (20%) creo que no habrá problemas. 
                                span=0.3)
    log_message("Selección de HVGs completada considerando efecto lote:\n")
except Exception as e2:
    log_message(f"Error capturado al intentar seleccionar HVGs corrigiendo por efecto lote: {e2}")
    sc.pp.highly_variable_genes(amp, layer="counts",  # El input son conteos
                                n_top_genes=top_genes,  # Top 20% HVGs
                                flavor="seurat_v3",  # Creo que es más adecuado para mi caso donde hay pacientes casos y controles. Para cada paciente, este método hace un ranking de x HVGs. Tras  ello, obtiene, para cada HVG, la mediana del ranking donde ha quedado. Finalmente, obtiene un ranking de HVGs donde se ordenan por mediana (de menor (más variable) a mayor (menos variable)), seleccionándose los x primeros HVGs. Este método tiene preferencia por HVGs que son consistentes, pudiendo ignorar HVG especificos de condición. Aún así, como el nº de HVGs es bastante alto (20%) creo que no habrá problemas. 
                                span=0.3)
    log_message("Selección de HVGs completada sin considerar efecto lote:\n")
print(amp)


# ============================================
# PCA
# ============================================
log_message("Obteniendo PCA...")
sc.pp.pca(amp, n_comps=50)
log_message("PCA completado :\n")
print(amp)


# ============================================
# CHECKPOINT 1
# ============================================
log_message(f"CHECKPOINT 1: Guardando el amp procesado en {output_amp}...\n")
print(amp)
amp.write_h5ad(output_amp, compression="gzip")
log_message("Listo!")


# ============================================
# CLUSTERING
# ============================================
# Calculo de vecinos
log_message("Obteniendo vecinos más cercanos y gráfico de vecindades...")
sc.pp.neighbors(amp, n_neighbors=30, n_pcs=50, metric="manhattan")
log_message("Cálculo de vecindades completado:\n")
print(amp)

# Calculo de clústeres
log_message("Ejecutando Leiden clustering a res = 0.1...")
sc.tl.leiden(amp, resolution=0.1, flavor="leidenalg", key_added="leiden_0_10")
log_message("Ejecutando Leiden clustering a res = 0.25...")
sc.tl.leiden(amp, resolution=0.25, flavor="leidenalg", key_added="leiden_0_25")
log_message("Ejecutando Leiden clustering a res = 1...")
sc.tl.leiden(amp, resolution=1, flavor="leidenalg", key_added="leiden_1_00")
log_message("Clustering completado:\n")
print(amp)

# ============================================
# CHECKPOINT 2
# ============================================
log_message(f"CHECKPOINT 2: Guardando el amp procesado en {output_amp}\n...")
print(amp)
amp.write_h5ad(output_amp, compression="gzip")
log_message("Listo!")


# ============================================
# tSNE
# ============================================
log_message("Ejecutando tSNE...")
sc.tl.tsne(amp, n_pcs=50, perplexity=500, early_exaggeration=20)
log_message("tSNE completado:\n")
print(amp)


# ============================================
# GUARDAR EL OBJETO h5ad
# ============================================
log_message(f"Guardando el amp procesado en {output_amp}...")
amp.write_h5ad(output_amp, compression="gzip")
log_message("Listo!")
