#!/bin/bash
#SBATCH --job-name=labeltransfer_symph_donor_job                   # Nombre del trabajo
#SBATCH --output=08-server_outputs/labeltransfer_symph_donor_job%J.out # Archivo de salida
#SBATCH --error=08-server_outputs/labeltransfer_symph_donor_job%J.err  # Archivo de errores
#SBATCH --ntasks=1                             # Nº de CPUs (en cualquier nodo)
#SBATCH --cpus-per-task=8                      # Nº de CPUs (en un mismo nodo)
#SBATCH --mem=1800gb                             # Memoria RAM total que necesita el trabajo
#SBATCH --time=4-00:00:00                        # Duración
#SBATCH --constraint=cal                       # Tipo de nodos donde se ejecutará (cal = cualquiera menos GPU)

INPUT_AMP=$1 # Ejemplo. 05-data_objects/04-amp_corrected.h5ad
INPUT_REF=$2 # Ejemplo. 05-data_objects/03-ref_corrected.h5ad
OUTPUT_AMP=$3 # Ejemplo. 05-data_objects/05-amp_annotated.h5ad
SAVE_PATH=$4 # Ejemplo. 06-results/01-plots/

SCRIPT_FILE=02-scripts/06_01_LabelTransfer_Symph_donorlvl.py

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
python -u $SCRIPT_FILE $INPUT_AMP $INPUT_REF $OUTPUT_AMP $SAVE_PATH
if [ $? -ne 0 ]; then
    echo "Error al ejecutar el script de Python"
    exit 1
fi

echo "Trabajo completado exitosamente."
# Mostrar la hora de finalización
echo "Hora de finalización: $(date)"
