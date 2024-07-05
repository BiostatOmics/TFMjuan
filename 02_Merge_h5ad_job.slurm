#!/bin/bash
#SBATCH --job-name=merge_job                   # Nombre del trabajo
#SBATCH --output=08-server_outputs/merge%J.out # Archivo de salida
#SBATCH --error=08-server_outputs/merge%J.err  # Archivo de errores
#SBATCH --ntasks=1                             # Nº de CPUs (en cualquier nodo)
#SBATCH --cpus-per-task=2                      # Nº de CPUs (en un mismo nodo)
#SBATCH --mem=512gb                             # Memoria RAM total que necesita el trabajo
#SBATCH --time=04:00:00                        # Duración
#SBATCH --constraint=cal                       # Tipo de nodos donde se ejecutará (cal = cualquiera menos GPU)

INPUT_FILE1=$1 # Ejemplo (debería ser). 03-raw_data/neurons_sub.h5ad
INPUT_FILE2=$2 # Ejemplo (debería ser). 03-raw_data/non_neurons_sub.h5ad
OUTPUT_FILE=$3 # Ejemplo. 03-raw_data/reference_sub.h5ad
SCRIPT_FILE=02-scripts/03_Merge_h5ad.py

# Activar conda
echo "Activando Conda..."
module load miniconda/3_py10
if [ $? -ne 0 ]; then
    echo "Error al activar conda"
    exit 1
fi

# Ejecutar el script de Python
echo "Ejecutando el script de Python..."
python -u $SCRIPT_FILE $INPUT_FILE1 $INPUT_FILE2 $OUTPUT_FILE
if [ $? -ne 0 ]; then
    echo "Error al ejecutar el script de Python"
    exit 1
fi

echo "Trabajo completado exitosamente."
