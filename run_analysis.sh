#!/bin/bash
# run_analysis.sh
# Script to run all IBD stromal cell analysis scripts in sequence

# Set directory to script location
cd "$(dirname "$0")"

# Create necessary directories
mkdir -p data results figures

echo "=== IBD Stromal Cell Analysis ==="
echo "Running all analysis scripts in sequence..."

# Run preprocessing
echo -e "\n[1/5] Running preprocessing..."
Rscript scripts/01-preprocessing.R
if [ $? -ne 0 ]; then
  echo "Error in preprocessing step. Stopping."
  exit 1
fi

# Run stromal subtyping
echo -e "\n[2/5] Running stromal subtyping..."
Rscript scripts/02-stromal_subtyping.R
if [ $? -ne 0 ]; then
  echo "Error in stromal subtyping step. Stopping."
  exit 1
fi

# Run pathway correlation analysis
echo -e "\n[3/5] Running pathway correlation analysis..."
Rscript scripts/03-pathway_correlation_analysis.R
if [ $? -ne 0 ]; then
  echo "Error in pathway correlation analysis step. Stopping."
  exit 1
fi

# Run receptor-ligand analysis
echo -e "\n[4/5] Running receptor-ligand analysis..."
Rscript scripts/04-receptor_ligand_analysis.R
if [ $? -ne 0 ]; then
  echo "Error in receptor-ligand analysis step. Stopping."
  exit 1
fi

# Run statistical tests
echo -e "\n[5/5] Running statistical tests..."
Rscript scripts/05-statistical_tests.R
if [ $? -ne 0 ]; then
  echo "Error in statistical tests step. Stopping."
  exit 1
fi

echo -e "\nAll analyses completed successfully!"
echo "Results are available in the results/ directory"
echo "Figures are available in the figures/ directory" 