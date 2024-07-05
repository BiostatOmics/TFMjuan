#!/bin/bash
#SBATCH --job-name=preprocess_job                   # Nombre del trabajo
#SBATCH --output=08-server_outputs/preprocess%J.out # Archivo de salida
#SBATCH --error=08-server_outputs/preprocess%J.err  # Archivo de errores
#SBATCH --ntasks=1                             # Nº de CPUs (en cualquier nodo)
#SBATCH --cpus-per-task=8                      # Nº de CPUs (en un mismo nodo)
#SBATCH --mem=1850gb                             # Memoria RAM total que necesita el trabajo
#SBATCH --time=6-12:00:00                        # Duración
#SBATCH --constraint=cal                       # Tipo de nodos donde se ejecutará (cal = cualquiera menos GPU)

INPUT_FILE=$1 # Ejemplo. 03-raw_data/01-reference.h5ad
OUTPUT_FILE=$2 # Ejemplo. 05-data_objects/02-reference_processed.h5ad

SCRIPT_FILE=02-scripts/03_02_REF_Preprocess.py


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
python -u $SCRIPT_FILE $INPUT_FILE $OUTPUT_FILE
if [ $? -ne 0 ]; then
    echo "Error al ejecutar el script de Python"
    exit 1
fi

echo "Trabajo completado exitosamente."
# Mostrar la hora de finalización
echo "Hora de finalización: $(date)"
