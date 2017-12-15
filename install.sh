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
# conda env create -f envs/anaconda.nanopolish-env.yaml

# Install albacore
if [[ -z "$(which read_fast5_basecaller.py)" ]]
then
  echo "Installing Albacore"

  sudo apt-get update
  sudo apt-get install wget
  wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | sudo apt-key add -
  echo "deb http://mirror.oxfordnanoportal.com/apt trusty-stable non-free" | sudo tee /etc/apt/sources.list.d/nanoporetech.sources.list
  sudo apt-get update
fi

  # Edit this line with the appropriate file path to Albacore deb
if [[ "$(uname -s)" -eq "Linux" ]]
then

  if [[ -z "$(which read_fast5_basecaller.py)" ]]
  then
    echo "Installing Albacore"
    sudo dpkg -i ~/Downloads/python3-ont-albacore_2.0.2-1-xenial_all.deb
    sudo apt-get -f install
  else
    echo "albacore already installed"
  fi

else

  echo "Please install Albacore before running pipeline"

fi

# Install nanopolish
if [[ -z "$(which nanopolish)" ]]
then
  cd && git clone --recursive https://github.com/jts/nanopolish.git && cd nanopolish && make
  export $PATH=$PATH:~/nanopolish/
else
  echo "nanopolish already installed"
fi
