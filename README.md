Code and ideation for the [virtual cell competition](https://virtualcellchallenge.org/). 

[Background information](https://github.com/KunaalAgarwal/virtual_cell_background), [team](https://virtualcellchallenge.org/app). 


**Usage**:

Most of our implementations will begin as Jupyter Notebooks. Thus you should start by creating a virtual environment for the project (don't commit this): 

```python -m venv <environment_name>```

Activate it with ```. <environment_name>/bin/activate``` (Linux), ```source <environment_name>/bin/activate``` (MacOS)

Install packages: ```pip install requirements.txt```
- If you add any new packages please update the requirements.txt so everyone will be able to run your notebook (```pip freeze > requirements.txt```)

For large data files we will be using git large file storage (LFS). If you generate a csv of results or add in a new data file ensure you ```git lfs track <file_path>``` prior to committing. 

For any AI workflows (hooks, etc.) please commit these if you think they would be useful to the group. 