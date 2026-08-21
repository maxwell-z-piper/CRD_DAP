# Local installation
Clone the github repo to your local device:
```bash
cd /Path/to/place/CRD_DAP
git clone https://github.com/maxwell-z-piper/CRD_DAP.git
```

Then, from Terminal, make a conda environment using python 3.12:
```bash
conda create -n crd_dap python=3.12
conda activate crd_dap
```

Within this environment, cd into CRD_DAP. Upgrade pip and install all the necessary packages. Then, run the pytest command to ensure everything is installed correctly.

```bash
cd CRD_DAP
python -m pip install --upgrade pip
python -m pip install -e ".[science,dev]"

python -m pytest -q
```

If the tests pass successfully, you are ready to run the pipeline. Remember to duplicate the target_config_template.py script and populate it with your target's specific information, and to rename it. This is needed to run every step of the pipeline as 
```bash
cd scripts
python *specific pipeline script*.py --config */path/to/specific config*.py
```
