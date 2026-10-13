# Origin of `main.py`

`main.py` is the HGNN training script distributed with SigXTalk
(<https://github.com/LithiumHou/SigXTalk>, `vignettes/main.py`), which is released under the MIT License
(Copyright (c) 2025 Jiawen Hou, <https://github.com/LithiumHou/SigXTalk/blob/master/LICENSE>). Please cite the SigXTalk paper when using it.

The copy in this folder is the version that was used for the reported results. It sets the random seed for Python, NumPy and
PyTorch (`set_global_seed(args.seed)`, default seed 42) so that the HGNN step is reproducible.
