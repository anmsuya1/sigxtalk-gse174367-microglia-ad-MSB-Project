# Python environment for the hypergraph neural network (HGNN)

The pipeline starts `hgnn/main.py` from R (`SigXTalkR::Run_py_script`, which uses `reticulate` 1.47.0) inside a conda
environment called **`SigXTalk_py`**. The environment name is hard-coded in the pipeline (`conda_env = "SigXTalk_py"`).

## Required versions (from the SigXTalk instructions, <https://github.com/LithiumHou/SigXTalk>)

| Component | Version |
|---|---|
| Python | 3.8.x |
| PyTorch | 1.13.1 (CPU build; CUDA 11.7 build `1.13.1+cu117` exists for Linux/Windows) |
| pandas, numpy | imported by `main.py`; no version is prescribed by SigXTalk |
| SigXTalkpy | the Python module in `pythoncodes/` of the SigXTalk repository (newest commit at install time) |

## Set-up

```bash
conda create -n SigXTalk_py python=3.8        # or: mamba create -n SigXTalk_py python=3.8
conda activate SigXTalk_py

# PyTorch 1.13.1 - choose ONE line
pip install torch==1.13.1                                                                  # macOS (used here; CPU only)
pip install torch==1.13.1+cpu  --extra-index-url https://download.pytorch.org/whl/cpu      # Linux / Windows, CPU
pip install torch==1.13.1+cu117 --extra-index-url https://download.pytorch.org/whl/cu117   # Linux / Windows, NVIDIA GPU

pip install pandas numpy

# the HGNN module (option 1: straight from GitHub)
pip install git+https://github.com/LithiumHou/SigXTalk.git#subdirectory=pythoncodes
# option 2: clone the repository, then   cd SigXTalk/pythoncodes && pip install .
```

Check that it works:

```bash
python -c "import torch, pandas, numpy; from SigXTalkpy import predictor, training, preprocessing; print(torch.__version__)"
```

`main.py` imports `torch`, `pandas`, `numpy` and the `SigXTalkpy` package (`predictor`, `training`, `preprocessing`). On Apple
Silicon there is no CUDA, so the HGNN runs on the CPU (`--device cpu`), with random seed 42.

## What was not recorded

The versions of `pandas`, `numpy` and `SigXTalkpy` that were actually installed when the results were produced were not saved.
On the machine that ran the pipeline, record them and commit the files:

```bash
conda activate SigXTalk_py
python --version
conda env export -n SigXTalk_py > environment/sigxtalk_py_full.yml
pip freeze > environment/sigxtalk_py_requirements.txt
```
