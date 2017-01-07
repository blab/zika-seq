# Pipelines for Zika sequencing

#### Trevor Bedford<sup>1</sup>

<sup>1</sup>Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, WA, USA

## Data sync

Primary sequencing data lives on the Rhino FHCRC cluster at:

    /fh/fast/bedford_t/zika-seq/data/

And locally on the Meristem drive at:

    /Volumes/Meristem/data/

To sync Meristem to Rhino, run:

    rsync -azP <username>@rhino.fhcrc.org:/fh/fast/bedford_t/zika-seq/data/ /Volumes/Meristem/data/

Replacing `<username>` with your username. This `data/` directory is assumed to follow [a particular schema](../data/).

## Bioinformatic pipeline

Here, we use the ZiBRA project bioinformatic pipeline at [zibraproject/zika-pipeline](https://github.com/zibraproject/zika-pipeline/). This pipeline is instantiated in the Docker image [zibra/zibra](https://hub.docker.com/r/zibra/zibra/). Data processing is done entirely using Docker.

### Data volume

Create a named data volume that mirrors local `data/` to `data/` within container:

    docker create --name zibra-data -v /Volumes/Meristem/data:/data zibra/zibra    

This is to get data into the Docker container. Note that the path to local directory has to be an absolute path. [Notes on data schema are here](../data/).

### Samples volume

Create a named data volume that mirrors local sample metadata `samples/` to `samples/` within container:

    docker create --name zibra-samples -v /Volumes/Meristem/samples:/samples zibra/zibra

This is to get sample metadata into the Docker container. Note that the path to local directory has to be an absolute path. [Notes on metadata schema are here](../samples/).

### Build volume

Create a named data volume that mirrors local `build/` to `build/` within container:

    docker create --name zibra-build -v /Volumes/Meristem/build:/build zibra/zibra

This is to get data out of the Docker container. Note that the path to local directory has to be an absolute path.

### Start

Enter docker image:

    docker run -t -i --volumes-from zibra-data --volumes-from zibra-samples --volumes-from zibra-build zibra/zibra /bin/bash

Run script:

    python scripts/pipeline.py
