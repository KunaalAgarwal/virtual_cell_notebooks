Code and ideation for the [virtual cell challenge](https://virtualcellchallenge.org/). 

[Background information](https://github.com/KunaalAgarwal/virtual_cell_background), [team](https://virtualcellchallenge.org/app). 

# Setup

```git clone https://github.com/KunaalAgarwal/virtual_cell_notebooks.git```

**[uv](https://docs.astral.sh/uv/getting-started/features/)**: alternative python package manager (alternative to pip and the like)
- Singular tool which replaces pip, virtual environments, etc.

Install uv if haven't already: https://docs.astral.sh/uv/getting-started/installation/

Install dependencies with ```uv sync```

Add new packages ```uv add <package 1> <package 2> ...<package n>```
- After adding a new package simply version the uv.lock and pyproject.toml files in your next commit. 

**[vcc](https://vcc-cli-wiki.virtualcellchallenge.org/)**: competition specific package.

Login into the vcc command line with your account: ```vcc login --token-stdin``` 
- Insert the API key (https://virtualcellchallenge.org/app > Credentials > Generate API Key)

After logged in ```vcc whoami``` should return your information

For large data files  we will be using git large file storage (LFS). If you generate a csv of results or add in a new data file ensure it is tracked ```git lfs -ls-files``` and if needed ```git lfs track <file_path>``` prior to committing. 

For any AI workflows (hooks, etc.) please commit these if you think they would be useful to the group. 
- \vcc skill is including in the directory and can be used for by agents

# Contributing

Add your notebooks to the src/ directory and add a brief summary statement at the top so others can quickly get acclimated. Feel free to create a subdirectory in src as your name so you can keep your files siloed. 

Run with via vscode jupyter notebook interface via .venv/bin/python environment. 

Run with jupyterlab cli: 
- ```uv run juypter lab --no-browser src/```
- Grab the url from the result of this command and paste in browser