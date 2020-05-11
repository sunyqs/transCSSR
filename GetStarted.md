# transCSSR

**GET READY:** install the required Python dependencies (see below)

**GET STARTED:** you can start playing with the derivation of epsilon-machines through the Jupyter code 'demo_CSSR.ipynb'. It is set up to derive the epsilon-machine for the 'even process' (a binary process of 1s and 0s, with the restriction that it will always have an even number of consecutive 1s). It displays the visuals, does cross-validation, and calculates (1) the statistical complexity C_mu, (2) the entropy rate h_mu, (3) the excess entropy E (aka 'predictive information').


**AUTOMATE INPUT:** the scripts called 'run_cssr..._on_csv' allow to read in sequence data directly from csv files. The script 
*   `run_cssr4eM_on_csv` derives a separate epsilon-machine for each input column, and  
*   `run_cssr4eT_on_csv` derives an epsilon-transducer for each input-output pair. 

The repository transCSSR here contains a csv folder, including `Yt_test.csv` and `Xt_test.csv` with three sequences. By default,
*   `run_cssr4eM_on_csv` targets `csv/Yt_test.csv` and generates files at `output_cssr_yt_test_L1`, with the measures for the respective three epsilon-machines for word-length 1.
*   `run_cssr4eT_on_csv` targets `csv/full_Xt.csv` & `csv/full_Yt.csv` and generates files at `output_trans_L1`



**PYTHON DEPENDENCIES:** This implementation requires the following dependencies for Python 3.7+:
* [numpy](http://www.numpy.org)
* [scipy](http://www.scipy.org)
* [pandas](http://pandas.pydata.org)
* [igraph](http://igraph.org/python/)
* [pylab](http://wiki.scipy.org/PyLab)
* [matplotlib](http://matplotlib.org)
These packages may be installed individually. A useful collection of these and related libraries for scientific computing with Python may be installed using [Enthought](https://store.enthought.com) or [Anaconda](https://www.continuum.io/downloads).
Visualization of the epsilon-machines requires [Graphviz](http://graphviz.org) or similar software for reading [dot](http://en.wikipedia.org/wiki/DOT_(graph_description_language)) files. To use Graphviz within a Juptyer notebook using the Anaconda distribution of Python, install the `python-graphviz` package using:

```
conda install python-graphviz
```

**NOTE:** This version of transCSSR is for use with Python 3.7+. A legacy version for use with Python 2.7 is hosted [here](https://github.com/ddarmon/transCSSR2).
