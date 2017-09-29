#!/bin/bash
#
# Install zika-seq requirements
#
# Download and install miniconda
if [[ -z "$(which conda)" ]]
then
  echo "Installing miniconda"

  if [["$(uname -s)" -eq "Darwin" ]]
  then
    CONDA_SCRIPT=Miniconda3-latest-MacOSX-x86_64.sh
  else
    CONDA_SCRIPT=Miniconda3-latest-Linux-x86_64.sh
  fi

  CONDA_DIR=$HOME/miniconda3
  wget https://repo.continuum.io/miniconda/$CONDA_SCRIPT
  bash $CONDA_SCRIPT -b -p $CONDA_DIR

  CONDA_BIN_DIR=$CONDA_DIR/bin
  PATH=$CONDA_BIN_DIR:$PATH

  rm -rf $CONDA_SCRIPT
else
  echo "miniconda already installed"
fi

# Environment in which snakemake will be run
conda env create -f envs/anaconda.snakemake-env.yaml
conda env create -f envs/anaconda.pipeline-env.yaml
conda env create -f envs/anaconda.nanopolish-env.yaml
