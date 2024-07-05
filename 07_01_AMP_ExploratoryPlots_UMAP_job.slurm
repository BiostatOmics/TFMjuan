#!/bin/bash
#SBATCH --job-name=exploratory_umap_job                   # Nombre del trabajo
#SBATCH --output=08-server_outputs/exploratory_umap_job%J.out # Archivo de salida
#SBATCH --error=08-server_outputs/exploratory_umap_job%J.err  # Archivo de errores
#SBATCH --ntasks=1                             # Nº de CPUs (en cualquier nodo)
#SBATCH --cpus-per-task=2                      # Nº de CPUs (en un mismo nodo)
#SBATCH --mem=800gb                             # Memoria RAM total que necesita el trabajo
#SBATCH --time=4-00:00:00                        # Duración
#SBATCH --constraint=cal                       # Tipo de nodos donde se ejecutará (cal = cualquiera menos GPU)

INPUT_FILE=$1 # Ejemplo. 05-data_objects/01-amp_processed.h5ad
SAVE_PATH=$2 # Ejemplo. 06-results/01-plots/

SCRIPT_FILE=02-scripts/07_01_AMP_ExploratoryPlots_UMAP.py

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
python -u $SCRIPT_FILE $INPUT_FILE $SAVE_PATH
if [ $? -ne 0 ]; then
    echo "Error al ejecutar el script de Python"
    exit 1
fi

echo "Trabajo completado exitosamente."
# Mostrar la hora de finalización
echo "Hora de finalización: $(date)"
