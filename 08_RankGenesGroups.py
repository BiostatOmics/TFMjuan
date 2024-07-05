# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 20:27:19 2024

DE Analysis (Rank Genes Groups)

@author: Juan Tejedor
"""

# ============================================
# CARGA DE MÓDULOS
# ============================================
import argparse
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
parser = argparse.ArgumentParser(description="Plots")
parser.add_argument("input_amp", type=str)
parser.add_argument("output_query", type=str)
args = parser.parse_args()

input_amp = args.input_amp
output_query = args.output_query

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
convert = ["path_braak_lb", "path_braak_nft"]
for var in convert:
    amp.obs[var] = amp.obs[var].astype("category")

# ============================================
# ANÁLISIS DE DEGS: MARKER GENES
# ============================================
# Antes elimino las células que no son caso ni control para simplificar el DE
log_message("Eliminado células other...")
amp_sub = amp[amp.obs["case_control"].isin(["Case", "Control"])]

# Calculo DEGs
log_message("Calculando DEGs...")
conditions = ["leiden_0_10_corrected", "leiden_0_25_corrected",
              "leiden_1_00_corrected", "case_control",
              "path_braak_lb", "path_braak_nft",
              "supercluster_term"]
log_message("Condiciones a las que calcular el DE:\n")
print(conditions)

for condition in conditions:
    # Saltar la ejecución del análisis DE si alguna de las categorías tiene 1 sola observación
    if not (amp_sub.obs[condition].value_counts() <= 1).any():
        log_message(f"Calculando DEGs para {condition}...")
        sc.tl.rank_genes_groups(amp_sub, groupby=condition, reference="rest",
                                method="wilcoxon", pts=True,
                                key_added=f"rank_genes_groups_{condition}")

# ============================================
# GUARDAR EL OBJETO h5ad
# ============================================
log_message(f"Guardando el amp procesado en {output_query}...")
amp_sub.write_h5ad(output_query, compression="gzip")
log_message("Listo!")
