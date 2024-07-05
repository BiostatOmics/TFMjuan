# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 15:36:11 2024

Merging

@author: Juan Tejedor
"""

# ============================================
# CARGA DE MÓDULOS
# ============================================
import argparse
import anndata as ad

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.curdir, '02-scripts')))
from my_utils import log_message

# ============================================
# PARSING DE CLI A PYTHON
# ============================================

# Creo el objeto parseador
parser = argparse.ArgumentParser(description="Merging")

# Establezco argumentos a parsear
parser.add_argument("input_h5ad_1", type=str)
parser.add_argument("input_h5ad_2", type=str)
parser.add_argument("output_h5ad", type=str)

# Extraigo los argumentos del CLI a un objeto de Python
args = parser.parse_args()

# Recupero argumentos como variables
input_h5ad_1 = args.input_h5ad_1
input_h5ad_2 = args.input_h5ad_2
output_h5ad = args.output_h5ad

# ============================================
# MAIN
# ============================================

# Cargo ambos h5ad y veo sus atributos
log_message(f"Cargando h5ad_1 desde {input_h5ad_1}...")
h5ad_1 = ad.read_h5ad(input_h5ad_1)
h5ad_1
h5ad_1.X
h5ad_1.obs
h5ad_1.var
h5ad_1.uns
h5ad_1.obsm

log_message(f"Cargando h5ad_2 desde {input_h5ad_2}...")
h5ad_2 = ad.read_h5ad(input_h5ad_2)
h5ad_2
h5ad_2.X
h5ad_2.obs
h5ad_2.var
h5ad_2.uns
h5ad_2.obsm

# Fusiono ambos h5ad
log_message("Fusionando h5ad_1 y h5ad_2...")
reference = ad.concat([h5ad_1, h5ad_2],
                      label="dataset_title",
                      keys=["neuronal", "non_neuronal"],
                      join="outer",
                      merge="unique",
                      uns_merge="unique")

# Añado los objetos que quedaron sin fusionar a "uns"
log_message("Añadiendo metadatos adicionales...")
reference.uns['citation_neurons'] = h5ad_1.uns['citation']
reference.uns['title_neurons'] = h5ad_1.uns['title']

reference.uns['citation_non_neurons'] = h5ad_2.uns['citation']
reference.uns['title_non_neurons'] = h5ad_2.uns['title']

reference.uns['X_UMAP_neurons'] = h5ad_1.obsm['X_UMAP']
reference.uns['X_X_tSNE_neurons'] = h5ad_1.obsm['X_tSNE']

reference.uns['X_UMAP_non_neurons'] = h5ad_2.obsm['X_UMAP']
reference.uns['X_tSNE_non_neurons'] = h5ad_2.obsm['X_tSNE']

reference
reference.X
reference.obs
reference.var
reference.uns
reference.obsm

# Guardo el h5ad fusionado
log_message(f"Guardando el h5ad fusionado en {output_h5ad}")
reference.write_h5ad(output_h5ad, compression="gzip")
log_message("Listo!")
