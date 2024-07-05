#!/bin/bash
#SBATCH --job-name=manual_annotation_job                   # Nombre del trabajo
#SBATCH --output=08-server_outputs/manual_annotation_job%J.out # Archivo de salida
#SBATCH --error=08-server_outputs/manual_annotation_job%J.err  # Archivo de errores
#SBATCH --ntasks=1                             # Nº de CPUs (en cualquier nodo)
#SBATCH --cpus-per-task=4                      # Nº de CPUs (en un mismo nodo)
#SBATCH --mem=800gb                             # Memoria RAM total que necesita el trabajo
#SBATCH --time=4-00:00:00                        # Duración
#SBATCH --constraint=cal                       # Tipo de nodos donde se ejecutará (cal = cualquiera menos GPU)

INPUT_FILE=$1 # Ejemplo. 05-data_objects/05-amp_annotated_ingest_dissections.h5ad
OUTPUT_FILE=$2 # Ejemplo. 05-data_objects/06-amp_final_annotations.h5ad
RESULTS_PATH=$3 # Ejemplo. 06-results/02-results/annotations.tsv


SCRIPT_FILE=02-scripts/09_ManualAnnotation.py

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
python -u $SCRIPT_FILE $INPUT_FILE $OUTPUT_FILE $RESULTS_PATH
if [ $? -ne 0 ]; then
    echo "Error al ejecutar el script de Python"
    exit 1
fi

echo "Trabajo completado exitosamente."
# Mostrar la hora de finalización
echo "Hora de finalización: $(date)"
