Bootstrap: docker
From: condaforge/miniforge3:24.7.1-0

%labels
    Author Connor Murray
    Version v0.0.1
    Description "Singularity definition file for TOPCHeF Python analyses (pip-only install)"

%help
    This container provides a reproducible Python 3.11 environment
    for TOPCHeF analyses, installing dependencies via pip in the base environment.

%environment
    export PATH=/opt/conda/bin:$PATH
    export LC_ALL=C
    export LANG=C
    export PYTHONUNBUFFERED=1
    export MPLCONFIGDIR=/tmp/matplotlib

%post
    echo ">>> Installing system and Python dependencies for TOPCHeF analyses..."

    # Install Python packages via pip
    pip3 install --upgrade \
        numpy \
        pandas \
        scipy \
        seaborn \
        matplotlib \
        tqdm \
        statsmodels \
        scikit-learn \
        joblib \
        tensorqtl \
        pandas_plink \
        pydeseq2 \
        qtl

    # Clean up cache to reduce image size
    pip3 cache purge

    echo ">>> TOPCHeF Python environment setup complete."
