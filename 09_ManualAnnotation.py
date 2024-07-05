# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 12:29:22 2024

Manual Annotation

@author: Juan Andrés Tejedor Serrano
"""

# ============================================
# CARGA DE MÓDULOS
# ============================================
import argparse
import pandas as pd
import anndata as ad
import scanpy as sc
import re

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.curdir, '02-scripts')))
from my_utils import log_message

# ============================================
# PARSING DE CLI A PYTHON
# ============================================

# Creo el objeto parseador
parser = argparse.ArgumentParser(description="Plots")
parser.add_argument("input_amp", type=str)
parser.add_argument("output_amp", type=str)
parser.add_argument("output_tsv", type=str)

args = parser.parse_args()

input_amp = args.input_amp
output_amp = args.output_amp
output_tsv = args.output_tsv


# ============================================
# CONFIGURACIÓN INICIAL Y CARGA DE DATOS
# ============================================
sc.settings.verbosity = 4
sc.settings.n_jobs = 8


log_message(f"Cargando amp desde {input_amp}...")
amp = ad.read_h5ad(input_amp)
log_message(f"Datos cargados desde {input_amp}:\n")
print(amp)

# ============================================
# CREACIÓN DE OBJETO FINAL
# ============================================

# Convierto a categórica
convert = ["hoehn_yahr_stage", "path_braak_lb", "path_braak_nft"]
for var in convert:
    amp.obs[var] = amp.obs[var].astype("category")

# Eliminar "counts_per_million" (pesa mucho)
del amp.layers["cpm"]
log_message(f"Estado actual {amp}:\n")

# Crear un DataFrame a partir de las columnas de interés y convertir a enteros si es necesario
df = amp.obs[['leiden_0_10_corrected', 'leiden_0_25_corrected']].copy()
df = df.apply(pd.to_numeric)

# Definir los diccionarios de mapeo
mapping_0_10 = {
    0: 'Oligodendrocytes (OLGs)',
    2: 'Astrocytes (ASTs)',
    3: 'Microglia (MG)',
    4: 'Oligodendrocyte Precursor Cells (OPCs)',
    6: 'MGE-derived Interneurons (MGE-derived Interneurons)',
    7: 'Deep-layer Intratelencephalic PCP4+ Neurons (DL-ITC PCP4+ Neurons)',
    10: 'Deep-layer Corticothalamic and Layer 6b TLE4+ Neurons (DL-CTC-6b TLE4+ Neurons)',
    11: 'Endothelial Cells (ECs)',
    12: 'Deep-layer Near-projecting CRYM+ Neurons (DL-NPC CRYM+ Neurons)',
    13: 'Chandelier Cells (ChCs)',
    14: 'Upper-layer Intratelencephalic CUX2- Neurons (UL-ITC CUX2- Neurons)',
}

mapping_0_25 = {
    2: 'Upper-layer Intratelencephalic CUX2+ Neurons (UL-ITC CUX2+ Neurons)',
    7: 'Deep-layer Intratelencephalic PCP4- Neurons (DL-ITC PCP4- Neurons)',
    15: 'LAMP5+ Interneurons (LAMP5+ Interneurons)',
    10: 'CGE-derived Interneurons (CGE-derived Interneurons)',
    14: 'Neurovascular Support Cells (NVSCs)',
    17: 'Cerebrospinal Fluid-related Cells of the Dorsal Motor Nucleus of the Vagus (CSFCs-DMNX)',
    18: '**ELIMINAR**',
    16: 'DRD2+ Neurons (DRD2+ Neurons)',
}

# Para leiden_0_10_corrected
log_message("Asignando anotaciones a leiden_0_10_corrected...")
df_0_10 = df[['leiden_0_10_corrected']].copy()
df_0_10['annotations'] = df_0_10['leiden_0_10_corrected'].map(mapping_0_10)

# Para leiden_0_25_corrected
log_message("Asignando anotaciones a leiden_0_25_corrected...")
df_0_25 = df[['leiden_0_25_corrected']].copy()
df_0_25['annotations'] = df_0_25['leiden_0_25_corrected'].map(mapping_0_25)


# Combinar los DataFrames de anotaciones sin sobrescribir las ya existentes
df_combined = df_0_10.join(df_0_25, lsuffix='_0_10', rsuffix='_0_25')

# Imprimir las columnas para verificar que los nombres son correctos
print(df_combined.columns)

log_message("Resolviento anotaciones combinadas...")
# Resolver las anotaciones combinadas
def combine_annotations(row):
    for col in ['annotations_0_10', 'annotations_0_25']:
        if col in row and row[col] != 'Unknown' and pd.notna(row[col]):
            return row[col]
    return 'Unknown'

df_combined['final_annotations_desc'] = df_combined.apply(combine_annotations, axis=1)

log_message("Extrayendo abreviaturas de las anotaciones...")
# Extraer la abreviatura entre paréntesis
def extract_abbreviation(annotation):
    match = re.search(r'\(([^)]+)\)', annotation)
    if match:
        return match.group(1)
    return annotation

df_combined['annotations'] = df_combined['final_annotations_desc'].apply(extract_abbreviation)

# Verificar los resultados
print(df_combined.head())

log_message("Agregando columnas de anotaciones al objeto amp...")
# Agregar estas nuevas columnas de anotaciones a tu objeto AnnData
amp.obs['final_annotations_desc'] = df_combined['final_annotations_desc']
amp.obs['annotations'] = df_combined['annotations']

# Eliminar las células que tienen la anotación "**ELIMINAR**" del objeto AnnData
log_message("Eliminando células con anotación '**ELIMINAR**' del objeto AnnData...")
amp = amp[amp.obs['annotations'] != '**ELIMINAR**']


# ============================================
# GUARDAR EL OBJETO h5ad
# ============================================
log_message(f"Guardando el amp procesado en {output_amp}:\n{amp}")
amp.write_h5ad(output_amp, compression="gzip")

# Extraer anotaciones finales como .tsv
log_message(f"Guardando DataFrame final como archivo .tsv en {output_tsv}...")
df_final = amp.obs[['leiden_0_10_corrected', 'leiden_0_25_corrected', 'leiden_1_00_corrected', 'final_annotations_desc', 'annotations'. 'supercluster_term']]
df_final.to_csv(output_tsv, sep='\t', index_label="cell_id")

log_message("Listo!")
