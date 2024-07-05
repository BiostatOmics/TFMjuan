#!/bin/bash
#SBATCH --job-name=addcellmetadata_job                   # Nombre del trabajo
#SBATCH --output=08-server_outputs/addcellmetadata%J.out # Archivo de salida
#SBATCH --error=08-server_outputs/addcellmetadata%J.err  # Archivo de errores
#SBATCH --ntasks=1                             # Nº de CPUs (en cualquier nodo)
#SBATCH --cpus-per-task=2                      # Nº de CPUs (en un mismo nodo)
#SBATCH --mem=256gb                             # Memoria RAM total que necesita el trabajo
#SBATCH --time=1-00:00:00                        # Duración
#SBATCH --constraint=cal                       # Tipo de nodos donde se ejecutará (cal = cualquiera menos GPU)

INPUT_FILE=$1 # Ejemplo. 05-data_objects/01-amp_processed.h5ad
CELL_METADATA=$2 # Ejemplo. 04-meta_data/corrected_full_cell_metadata.tsv
OUTPUT_FILE=$3 # Ejemplo. 05-data_objects/01-amppd_processed.h5ad

SCRIPT_FILE=02-scripts/04_AMP_AddCellMetadata.py


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
python -u $SCRIPT_FILE $INPUT_FILE $CELL_METADATA $OUTPUT_FILE
if [ $? -ne 0 ]; then
    echo "Error al ejecutar el script de Python"
    exit 1
fi

echo "Trabajo completado exitosamente."
# Mostrar la hora de finalización
echo "Hora de finalización: $(date)"
