#!/bin/bash
#SBATCH --job-name=ingest_dissections_labeltransfer_job                   # Nombre del trabajo
#SBATCH --output=08-server_outputs/ingest_dissections_labeltransfer_job%J.out # Archivo de salida
#SBATCH --error=08-server_outputs/ingest_dissections_labeltransfer_job%J.err  # Archivo de errores
#SBATCH --ntasks=1                             # Nº de CPUs (en cualquier nodo)
#SBATCH --cpus-per-task=12                      # Nº de CPUs (en un mismo nodo)
#SBATCH --mem=1800gb                             # Memoria RAM total que necesita el trabajo
#SBATCH --time=7-00:00:00                        # Duración
#SBATCH --constraint=cal                       # Tipo de nodos donde se ejecutará (cal = cualquiera menos GPU)

INPUT_QUERY=$1 # Ejemplo. 05-data_objects/04-amppd_corrected.h5ad
INPUT_REF=$2 # Ejemplo. 01-reference.h5ad
OUTPUT_QUERY=$3 # Ejemplo. 05-data_objects/05-amp_annotated_ingest_dissections.h5ad
OUTPUT_REF=$4 # Ejemplo. 05-data_objects/04-ref_from_ingest_dissections.h5ad
SAVE_PATH=$5 # Ejemplo. 06-results/01-plots/exploratory_umaps_ingest_dissections

SCRIPT_FILE=02-scripts/06_02_LabelTransfer_Ingest_dissections.py

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
python -u $SCRIPT_FILE $INPUT_QUERY $INPUT_REF $OUTPUT_QUERY $OUTPUT_REF $SAVE_PATH
if [ $? -ne 0 ]; then
    echo "Error al ejecutar el script de Python"
    exit 1
fi

echo "Trabajo completado exitosamente."
# Mostrar la hora de finalización
echo "Hora de finalización: $(date)"
