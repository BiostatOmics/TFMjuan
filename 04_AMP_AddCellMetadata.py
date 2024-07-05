# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 16:23:38 2024

Cell Metadata Transfer

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

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.curdir, '02-scripts')))
from my_utils import log_message
# ============================================
# PARSING DE CLI A PYTHON
# ============================================

# Creo el objeto parseador
parser = argparse.ArgumentParser(description="Cell Metadata Transfer")

# Establezco argumentos a parsear
parser.add_argument("input_amp", type=str)
parser.add_argument("cell_metadata_path", type=str)
parser.add_argument("output_amp", type=str)

# Extraigo los argumentos del CLI a un objeto de Python
args = parser.parse_args()

# Recupero argumentos como variables
input_amp = args.input_amp
cell_metadata_path = args.cell_metadata_path
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

log_message(f"Cargando metadatos celulares desde {cell_metadata_path}...")
cell_metadata = pd.read_csv(cell_metadata_path, sep=",", header=0, index_col=0)
log_message(f"Datos cargados desde {cell_metadata_path}:\n")
print(cell_metadata)

# ============================================
# TRANSFERENCIA DE METADATOS CELULARES
# ============================================

log_message("Transfiriendo metadatos celulares...")
# Estos metadatos contienen los metadatos celulares detallados, las métricas de QC obtenidas en Seurat y la anotación celular obtenida por ingest.
# Identifico las columnas no presentes en amp
columns_to_transfer = cell_metadata.columns.difference(amp.obs.columns)
# Filtro la tabla con los metadatos para que sólo contenga estas columnas
cell_metadata_to_transfer = cell_metadata[columns_to_transfer]
# Filtro las células para que coincidan en ambos datasets
cell_metadata_to_transfer = cell_metadata_to_transfer[cell_metadata_to_transfer.index.isin(amp.obs.index)]
# Fusiono ambos metadatos transfiriendo sólo las columnas nuevas
amp.obs = pd.concat([amp.obs, cell_metadata_to_transfer], axis=1)
# Paso a categóricas
# Especifico las columnas a convertir
amp.obs["hoehn_yahr_stage"] = amp.obs["hoehn_yahr_stage"].astype('category')
amp.obs["path_braak_lb"] = amp.obs["path_braak_lb"].astype('category')
amp.obs["path_braak_nft"] = amp.obs["path_braak_nft"].astype('category')
amp.obs["barcodekey"] = amp.obs["barcodekey"].astype('category')


log_message("Metadatos celulares transferidos:\n")
print(amp)


# ============================================
# GUARDAR EL OBJETO h5ad
# ============================================
log_message(f"Guardando el amp procesado en {output_amp}...")
amp.write_h5ad(output_amp, compression="gzip")
log_message("Listo!")
