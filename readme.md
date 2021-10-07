# R Plotting Functions

## Overview

* A collection of R functions. New plotting functions will be added along with examples as I finish them
* Partially to help scientists and others who want to avoid proprietary software like Excel and GraphPad's [Prism](http://www.graphpad.com/scientific-software/prism/). This allows analysis software/code to be portable across labs and countries. 
	* Also encourages storage of data portable, plain text file formats (e.g. .csv) rather than binary formats (e.g. .xls).

## Features

* Commented R plotting functions (see `\functions`).
* Example scripts for each function (see `\examples`).
* Example output plots (using given functions) to compare against.
* Tested in R 2.15.2.

## Installation

### Git

Go to your personal LaTeX projects directory (or wherever) and clone the repository using the command below:

    git clone https://github.com/gaara42/R_plot_functions

### Manual

* Download the files using the GitHub .zip download option
* Unzip the files and rename the folder to whatever name you want, references are relative only below the root folder.

### Dependencies

* Download [R](http://www.r-project.org/)). 
	* Make sure to have all packages updated.
	* Several packages need to be added to get function to work

## Instructions

###Overview

* Tested in R 2.15.2, can check old version if people find compatibility issues.
* `\functions` contains plotting functions, with naming convention `view.PLOT_TYPE.R`.
* `\examples` are example plotting functions, with naming convention `analysis.PLOT_TYPE.EXAMPLE_NAME.R`.
* `\plots` are example plots, with naming convention `plot.PLOT_TYPE.EXAMPLE_NAME.jpg`.
	* JPEGs are used in the repository since they are widely supported, use `postscript()` function to get vector output, which is higher quality, especially if combining R with LaTeX.

###PCA Function

* Use to analyze datasets with many dimensions for a few conditions. For example, looking at [microarray](http://en.wikipedia.org/wiki/DNA_microarray) (gene expression) data for patients with or without cervical cancer.
* `\functions\view.pca.R` 
	* The main function `pcaplot` accepts both data.frame and string reference to a file. 
	* If given a file path, it will attempt to parse into a data.frame. 
		* Make sure file has headers and row names. If not, set `rownames` and `headers` to `FALSE` when calling `pcaplots`
* `\examples\analysis.pca.USA.R` is an example implementation usng the `USArrests` data.frame that is auto-loaded in 2.15.2
	* Note, `USArrests` is the correct input format for this function. 

## Changelog

### What's New

* Cleaned up PCA function and added description of inputs (2013.01.10)
* Added to GitHub (2013.01.10)

### In Development

* Clean-up code and comments, add additional R features.
* Codes for: box, bar, scatter, line, density and other types of common plots.
* More examples!

## License

Biafra Ahanonu

This program is free software: you can redistribute it and/or modify it under the terms of the [GNU](http://www.gnu.org/licenses/gpl.html). Attribution is appreciated, but not required, if parts of the software are used elsewhere.
