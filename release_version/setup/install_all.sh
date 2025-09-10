
#!/bin/bash
# Create a new environment with the given name (or default metasrna_env)
# Usage: ./install.sh [env_name]
# Example: ./install_all.sh meta_test_all
#          ./install_all.sh        # will use metasrna_all

set -euo pipefail

ENV_NAME="${1:-metasrna_all}"

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
$PKG create -n "$ENV_NAME" -y python=3.8 blast=2.14.0 mirdeep2 viennarna=1.8.5 seqtk=1.2 entrez-direct=16.2 matplotlib tensorflow=2.13.*

echo "[info] activating $ENV_NAME"
# source ~/.bashrc 2>/dev/null || true
source $(conda info --base)/etc/profile.d/conda.sh

conda activate "$ENV_NAME" || { echo "Error: failed to activate $ENV_NAME"; exit 1; }
export PATH="$CONDA_PREFIX/bin:$PATH"

echo "[info] installing dependencies"

$PKG install -y -c conda-forge numpy pandas seaborn pillow scikit-learn scipy brokenaxes intervaltree openpyxl
$PKG install -y -c conda-forge biopython scikit-bio pybedtools
$PKG install -y -c conda-forge ipython jupyter
$PKG install -y -c bioconda linearfold
# Levenshtein (string distance)
$PKG install -y -c conda-forge python-Levenshtein

# Optional: TensorFlow for neural net-based models
# $PKG install -y -c conda-forge "tensorflow>=2.15.0"
echo "[done] environment '$ENV_NAME' ready at: $CONDA_PREFIX"
echo "       activate with: conda activate $ENV_NAME"

###########################################################
# Configure ViennaRNA 1.8.4(RNAfold) and PERL5LIB fix
###########################################################
script_path="$(readlink -f "$0")"
dir="$(dirname "$script_path")"

target_path_vienna="$dir/ViennaRNA-1.8.4/install_dir/bin"

mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
mkdir -p "$CONDA_PREFIX/etc/conda/deactivate.d"

# activation script
cat > "$CONDA_PREFIX/etc/conda/activate.d/env_vars.sh" <<EOF
#!/bin/bash
# Save PATH only if not already saved
if [[ -z "\${OLD_PATH-}" ]]; then
    export OLD_PATH="\$PATH"
fi
# Prepend ViennaRNA bins
export PATH=$target_path_vienna:\$PATH

# Force Perl to use this env's libraries
export PERL5LIB=\$CONDA_PREFIX/lib/perl5/site_perl:\$PERL5LIB
EOF

# deactivation script
cat > "$CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh" <<EOF
#!/bin/bash
# Restore PATH only if OLD_PATH is set
if [[ -n "\${OLD_PATH-}" ]]; then
    export PATH="\$OLD_PATH"
    unset OLD_PATH
fi
# Reset PERL5LIB
unset PERL5LIB
EOF

chmod +x "$CONDA_PREFIX/etc/conda/activate.d/env_vars.sh"
chmod +x "$CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh"

echo "[info] Finished setup."
conda deactivate