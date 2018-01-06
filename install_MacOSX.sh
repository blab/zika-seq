#!/bin/bash
# Install zika-seq requirements -- MacOSX

if [[ -z "$(which wget)" ]]
then
  echo "Installing wget"
  brew install wget
  if [[ -z "$(which wget)" ]]
  then
    echo "Successfully intalled wget"
  else
    echo "Error installing wget: skipping for now."
    echo "Make sure that wget is properly installed if Albacore, Porechop, or Nanopolish fail to install."
  fi
else
  echo "Wget already installed"
fi

## Install Miniconda3
if [[ -z "$(which conda)" ]]
then
  echo "Installing Miniconda"
  CONDA_SCRIPT=Miniconda3-latest-MacOSX-x86_64.sh

  CONDA_DIR=$HOME/miniconda3
  wget https://repo.continuum.io/miniconda/$CONDA_SCRIPT
  bash $CONDA_SCRIPT -b -p $CONDA_DIR

  CONDA_BIN_DIR=$CONDA_DIR/bin
  export PATH=$CONDA_BIN_DIR:$PATH
  if [[ -z "$(which conda)" ]]
  then
    echo "Successfully installed Miniconda"
  else
    echo "Error installing Miniconda: skipping for now"
  fi
else
  echo "Miniconda already installed"
fi

## Create environment in which snakemake will be run
echo "Installing conda environment: zika-seq"
conda env create -f envs/anaconda.snakemake-env_MacOSX.yaml
echo "Installing conda environment: zika-seq_pipeline"
conda env create -f envs/anaconda.pipeline-env.yaml

## Install albacore
if [[ -z "$(which read_fast5_basecaller.py)" ]]
then
  echo "Installing Albacore"
  # Edit this line with path to Albacore wheel file downloaded from https://community.nanoporetech.com/downloads
  # We recommend creating a zika-seq/albacore directory and downloading to that directory
  pip3 install albacore/ont_albacore-xxx_x86_64.whl || echo "Please donwload Albacore wheel file from https://community.nanoporetech.com/downloads and update install_macOSX.sh line 38 as appropriate."
  if [[ -z "$(which read_fast5_basecaller.py)" ]]
  then
    echo "Successfully installed Albacore"
  else
    echo "Error installing Albacore: skipping for now"
  fi
else
  echo "Albacore already installed"
fi

## Install Porechop
if [[ -z "$(which porechop)" ]]
then
  git clone https://github.com/rrwick/Porechop.git && cd Porechop && python3 setup.py install
  if [[ -z "$(which porechop)" ]]
  then
    echo "Successfully installed Porechop"
  else
    echo "Error installing Porechop: skipping for now"
  fi
else
  echo "Porechop already installed"
fi

## Install Nanopolish
if [[ -z "$(which nanopolish)" ]]
then
  git clone --recursive https://github.com/jts/nanopolish.git && cd nanopolish && make CXX=g++-7 CC=gcc-7 && cd ..
  NANOPOLISH_BIN_DIR=$PATH/nanopolish/bin
  export PATH=$NANOPOLISH_BIN_DIR:$PATH
  if [[ -z "$(which nanopolish)" ]]
  then
    echo "Successfully installed Nanopolish"
  else
    echo "Error installing Nanopolish: skipping for now"
  fi
else
  echo "Nanopolish already installed"
fi

rm -rf Porechop/
rm -rf $CONDA_SCRIPT
