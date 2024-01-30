# switch
This repo contains material on network switches and biological oscillations

This repository contains the reproducible code for the generation of data, analysis and figures in the manuscript "Mechanisms generating network switches and their role in circadian clocks". <!--([preprint](https://www.biorxiv.org/content/10.1101/2023.02.12.528191v1)).-->

To execute this code:

1. Clone this repository at a suitable location. This will place the code within a directory named **switch**
<!--2. Download all the simulated data from [here](https://www.zenodo.org/) (under the `results` folder) (or alternatively generate all the simulated data using the *main.py* script)-->
2. Make sure you have Python 3.7+ and the following libraries installed: `numpy`, `scipy`, `matplotlib`, `ddeint`
3. To reproduce the boxes and figures within, the Python files can now be executed *box1.py*, *box2.py,* ... within the project 

The Python script used to generate results and figures rely on objects stored in the `utils/` directory. The results for the bifurcation analyses from Box 6 have been generated with XPP-AUTO, and are stored in the `results/` directory. 

<!--To reproduce the figures, the Python files can now be executed in order
