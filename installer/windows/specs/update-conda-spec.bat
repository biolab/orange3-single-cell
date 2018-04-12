conda config --add channels conda-forge
conda remove -y -n orange --all
conda create -y -n orange python=3.6.* Orange3 keyring=9.0 scipy=0.18.1 numpy=1.12.1 orange3-bioinformatics
conda list -n orange --export --explicit --md5 > conda-spec.txt
