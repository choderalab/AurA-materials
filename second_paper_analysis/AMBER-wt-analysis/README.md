## Manifest
* `analysis_scripts` - all python code used to calculate specific distance distributions used in the data analysis
* `distance_analysis.ipynb` - [Jupyter Notebook](http://jupyter.org/) used to generate figures 

## AurA simulation data analysis 

This directory contains all of the python code to analysis the Folding@Home Simulation data, primarily relying on mdtraj and numpy. The ipython notebook used the output numpy array files to create the figures used in the publication, relying on matplotlib and seaborn. 

We suggested you install miniconda and then install dependencies with: 

```
conda install --yes jupyter seaborn scipy numpy matplotlib mdtraj
```

We will make the Folding@Home data available from the [Open Science Framework](https://osf.io/afg8h/) shortly. You can use the provided scripts to explore the data on your own.  
