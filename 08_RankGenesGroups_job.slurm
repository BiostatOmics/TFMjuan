#!/bin/bash
#SBATCH --job-name=rank_genes_job                   # Nombre del trabajo
#SBATCH --output=08-server_outputs/rank_genes_job%J.out # Archivo de salida
#SBATCH --error=08-server_outputs/rank_genes_job%J.err  # Archivo de errores
#SBATCH --ntasks=1                             # Nº de CPUs (en cualquier nodo)
#SBATCH --cpus-per-task=4                      # Nº de CPUs (en un mismo nodo)
#SBATCH --mem=1500gb                             # Memoria RAM total que necesita el trabajo
#SBATCH --time=6-00:00:00                        # Duración
#SBATCH --constraint=cal                       # Tipo de nodos donde se ejecutará (cal = cualquiera menos GPU)

INPUT_QUERY=$1 # Ejemplo. 05-data_objects/04-amppd_corrected.h5ad
OUTPUT_QUERY=$2 # Ejemplo. 05-data_objects/05-amp_annotated_ingest.h5ad

SCRIPT_FILE=02-scripts/08_RankGenesGroups.py

# Mostrar la hora de inicio
echo "Hora de inicio: $(date)"

# Activar conda
echo "Activando Conda..."
module load miniconda/3_py10
if [ $? -ne 0 ]; then
    echo "Error al activar conda"
    exit 1
fi

# Ejecutar el script de Python
echo "Ejecutando el script de Python..."
python -u $SCRIPT_FILE $INPUT_QUERY $OUTPUT_QUERY
if [ $? -ne 0 ]; then
    echo "Error al ejecutar el script de Python"
    exit 1
fi

echo "Trabajo completado exitosamente."
# Mostrar la hora de finalización
echo "Hora de finalización: $(date)"
