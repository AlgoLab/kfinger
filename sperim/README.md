# Experiments reproduction

All the experiments can be reproduced using the [`Snakefile`](./Snakefile) that is provided in this directory.

All the dependencies are managed using either `conda` or (preferably) `mamba`.
Please either install [miniconda](https://docs.conda.io/en/latest/miniconda.html) or
[mambaforge](https://github.com/conda-forge/miniforge#mambaforge) according to the instructions for your OS
(if it is not already installed).

Then, install all the dependencies in a new environment (here named `kfinger_sperim`) using the command:

``` shell
mamba env create -n kfinger_sperim -f environment.yml
conda activate kfinger_sperim
```

Replace `mamba` with `conda` in the first command if you are using miniconda.

Finally, you are able to reproduce the experiments with the command:

``` shell
snakemake -j4 -p --use-conda
```

where `4` is the number of cores the system might use  for the experiments (please increase/decrease it depending on your system).

