# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 11:42:44 2024

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
parser.add_argument("input_ref", type=str)
parser.add_argument("output_ref", type=str)

# Extraigo los argumentos del CLI a un objeto de Python
args = parser.parse_args()

# Recupero argumentos como variables
input_ref = args.input_ref
output_ref = args.output_ref

# ============================================
# CONFIGURACIÓN INICIAL Y CARGA DE DATOS
# ============================================

# Establezco ajustes de Scanpy
sc.settings.verbosity = 4
sc.settings.n_jobs = 8

# Cargo datos y veo sus atributos
log_message(f"Cargando ref desde {input_ref}...")
ref = ad.read_h5ad(input_ref)
log_message(f"Datos cargados desde {input_ref}:\n")
print(ref)


# ============================================
# CONTROL DE CALIDAD
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


# ============================================
# DETECCIÓN DE DOUBLETS
# ============================================

log_message("Detectando doublets usando scrublet...")
# Intentar corregir efecto lote por participante, de lo contrario, no corregir
try:
    ref = sc.pp.scrublet(ref, n_prin_comps=50, batch_key="sample_id",
                         copy=True, verbose=True)
    log_message("Detección completada corregiendo efecto lote por muestra")
except Exception as e1:
    log_message(f"Error capturado al intentar corregir efecto lote en scrublet por muestra: {e1}")
    ref = sc.pp.scrublet(ref, n_prin_comps=50,
                         copy=True, verbose=True)
    log_message("Detección completada sin corregir efecto lote")

log_message("Detección de doublets completada:\n")
print(ref)

# ============================================
# NORMALIZACIÓN
# ============================================
log_message("Matriz de expresión pre-normalización:\n")
print(ref.X.toarray())

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


# ============================================
# SELECCIÓN DE HVGs
# ============================================

# Selección de HVGs
log_message("Seleccionando HVGs...")
top_genes = np.floor(0.2*ref.shape[1]).astype(int)
try:
    sc.pp.highly_variable_genes(ref, layer="counts",  # El input son conteos
                                n_top_genes=top_genes,  # Top 20% HVGs
                                batch_key="sample_id",  # Considero efecto lote
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


# ============================================
# PCA
# ============================================
log_message("Obteniendo PCA...")
sc.pp.pca(ref, n_comps=50)
log_message("PCA completado :\n")
print(ref)


# ============================================
# CHECKPOINT 1
# ============================================
log_message(f"CHECKPOINT 1: Guardando el ref procesado en {output_ref}...\n")
print(ref)
ref.write_h5ad(output_ref, compression="gzip")
log_message("Listo!")


# ============================================
# CLUSTERING
# ============================================
# Calculo de vecinos
log_message("Obteniendo vecinos más cercanos y gráfico de vecindades...")
sc.pp.neighbors(ref, n_neighbors=30, n_pcs=50, metric="manhattan")
log_message("Cálculo de vecindades completado:\n")
print(ref)

# Calculo de clústeres
log_message("Ejecutando Leiden clustering a res = 0.1...")
sc.tl.leiden(ref, resolution=0.1, flavor="leidenalg", key_added="leiden_0_10")
log_message("Clustering completado:\n")
print(ref)

# ============================================
# CHECKPOINT 2
# ============================================
log_message(f"CHECKPOINT 2: Guardando el ref procesado en {output_ref}\n...")
print(ref)
ref.write_h5ad(output_ref, compression="gzip")
log_message("Listo!")


# ============================================
# tSNE
# ============================================
log_message("Ejecutando tSNE...")
sc.tl.tsne(ref, n_pcs=50, perplexity=500, early_exaggeration=20)
log_message("tSNE completado:\n")
print(ref)


# ============================================
# GUARDAR EL OBJETO h5ad
# ============================================
log_message(f"Guardando el ref procesado en {output_ref}...")
ref.write_h5ad(output_ref, compression="gzip")
log_message("Listo!")
