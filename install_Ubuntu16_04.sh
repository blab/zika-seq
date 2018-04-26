#!/bin/bash
#
# Install zika-seq requirements
#
# Download and install miniconda
if [[ -z "$(which conda)" ]]
then
  echo "Installing miniconda"

  CONDA_SCRIPT=Miniconda3-latest-Linux-x86_64.sh

  CONDA_DIR=$HOME/miniconda3
  wget https://repo.continuum.io/miniconda/$CONDA_SCRIPT
  bash $CONDA_SCRIPT -b -p $CONDA_DIR

  CONDA_BIN_DIR=$CONDA_DIR/bin
  export PATH=$CONDA_BIN_DIR:$PATH

  rm -rf $CONDA_SCRIPT
  if [[ -z "$(which conda)" ]]
  then
    echo "Successfully installed Miniconda"
  else
    echo "Error installing Miniconda: skipping for now"
  fi
else
  echo "miniconda already installed"
fi

# Environment in which snakemake will be run
echo "Installing conda environment: zika-seq"
conda env create -f envs/anaconda.snakemake-env.yaml
echo "Installing conda environment: zika-seq_pipeline"
conda env create -f envs/anaconda.pipeline-env.yaml

# Install Albacore
if [[ -z "$(which read_fast5_basecaller.py)" ]]
then
  echo "Installing Albacore"

  sudo apt-get update
  sudo apt-get install wget
  wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | sudo apt-key add -
  echo "deb http://mirror.oxfordnanoportal.com/apt trusty-stable non-free" | sudo tee /etc/apt/sources.list.d/nanoporetech.sources.list
  sudo apt-get update
  # Edit this line with the appropriate file path to Albacore deb
  sudo dpkg -i ~/Downloads/python3-ont-albacore_2.0.2-1-xenial_all.deb
  sudo apt-get -f install
else
  echo "Albacore already installed"
fi

# Install Nanopolish
if [[ -z "$(which nanopolish)" ]]
then
  cd && git clone --recursive https://github.com/jts/nanopolish.git && cd nanopolish && make
  export $PATH=$PATH:~/nanopolish/
  if [[ -z "$(which nanopolish)" ]]
  then
    echo "Successfully installed Nanopolish"
  else
    echo "Error installing Nanopolish: skipping for now"
  fi
else
  echo "Nanopolish already installed"
fi

# echo "Building data directories and downloading test data"
# mkdir -p data/usvi-library8-1d-2017-03-31/process/demux
# mkdir data/usvi-library8-1d-2017-03-31/raw_reads
# mkdir data/usvi-library8-1d-2017-03-31/basecalled_reads
# mkdir build
#
# # extract test data
# echo "Downloading and extracting example dataset"
# cd data/usvi-library8-1d-2017-03-31/raw_reads/ && wget https://s3.amazonaws.com/trvrb/zika_seq_example.tar.gz
# tar xvzf data/usvi-library8-1d-2017-03-31/raw_reads/zika_seq_example.tar.gz
# rm data/usvi-library8-1d-2017-03-31/raw_reads/zika_seq_example.tar.gz

echo "Done!"
