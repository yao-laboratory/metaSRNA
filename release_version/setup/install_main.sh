
#!/bin/bash
# Create a new environment with the given name (or default metasrna_env)
# Usage: ./install.sh [env_name]
# Example: ./install_main.sh meta_test_main
#          ./install_main.sh        # will use metasrna_main

set -euo pipefail

ENV_NAME="${1:-metasrna_main}"

# check for conda
if ! command -v conda >/dev/null 2>&1; then
  echo "Error: conda not found in PATH. Install Miniconda/Anaconda first."
  exit 1
fi

# initialize shell functions if needed
if ! type conda >/dev/null 2>&1; then
  conda init bash >/dev/null 2>&1 || true
  source ~/.bashrc 2>/dev/null || true
fi

# try to get mamba; fallback to conda
PKG=conda
if ! command -v mamba >/dev/null 2>&1; then
  echo "[info] mamba not found; trying to install in base..."
  conda install -n base -c conda-forge -y mamba || true
  hash -r
fi
if command -v mamba >/dev/null 2>&1; then
  PKG=mamba
fi
echo "[info] using package tool: $PKG"

# check if environment already exists
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "[warn] environment '$ENV_NAME' already exists. Deleting it first..."
  conda env remove -n "$ENV_NAME" -y
fi

echo "[info] creating environment: $ENV_NAME (python)"
$PKG create -n "$ENV_NAME" -y python=3.10 blast=2.14.0 seqtk=1.2 mirdeep2 linearfold

echo "[info] activating $ENV_NAME"
# source ~/.bashrc 2>/dev/null || true
source $(conda info --base)/etc/profile.d/conda.sh
conda activate "$ENV_NAME" || { echo "Error: failed to activate $ENV_NAME"; exit 1; }

echo "[info] installing dependencies"

$PKG install -y -c conda-forge numpy pandas openpyxl
$PKG install -y -c conda-forge biopython pybedtools
$PKG install -y -c conda-forge ipython jupyter
# $PKG install -y -c bioconda linearfold


echo "[done] environment '$ENV_NAME' ready at: $CONDA_PREFIX"
echo "       activate with: conda activate $ENV_NAME"

# Absolute path of the current script
script_path="$(readlink -f "$0")"
# Parent folder (folder)
dir="$(dirname "$script_path")"
# The bin folder under that parent (folder/bin)
target_path="$dir/viennarna_1.8.4_bin"
# ensure ~/.bashrc won't crash on nounset with /etc/bashrc
if ! grep -Fq ': "${BASHRCSOURCED:=1}"' ~/.bashrc; then
  sed -i '1i : "${BASHRCSOURCED:=1}"' ~/.bashrc
fi
# Add to PATH via ~/.bashrc if not already present
if ! grep -q "$target_path" ~/.bashrc; then
    echo "Adding $target_path to PATH in ~/.bashrc"
    echo "export PATH=\$PATH:$target_path" >> ~/.bashrc
else
    echo "$target_path already exists in ~/.bashrc"
fi
# Reload ~/.bashrc in this shell
echo "Reloading ~/.bashrc..."
# shellcheck source=/dev/null
source ~/.bashrc
conda deactivate

