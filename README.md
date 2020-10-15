## Description
The use of genetic variants as instrumental variables – an approach known 
as Mendelian randomization – is a popular epidemiological method for estimating 
the causal effect of an exposure (phenotype, biomarker, risk factor)
on a disease or health-related outcome from observational data. Instrumental 
variables must satisfy strong, often untestable assumptions, which means
that finding good genetic instruments among a large list of potential candidates 
is challenging. This difficulty is compounded by the fact that many
genetic variants influence more than one phenotype through different causal
pathways, a phenomenon called horizontal pleiotropy. This leads to errors
not only in estimating the magnitude of the causal effect, but also in inferring 
the direction of the putative causal link. We propose a Bayesian approach 
(BayesMR) that is a generalization of the Mendelian randomization technique in 
which we allow for pleiotropic effects and, crucially, for the possibility of 
reverse causation. The output of the method is a posterior distribution over the 
target causal effect, which provides an immediate and easily interpretable 
measure of the uncertainty in the estimation. More importantly, with the BayesMR
algorithm we can estimate the model evidence for both directions to determine 
how much more likely the expected direction is relative to the reverse direction.

## BayesMR

The data set contains source code implementing the BayesMR algorithm, which is 
described in the article titled "[Inferring the Direction of a Causal Link and 
Estimating Its Effect via a Bayesian Mendelian Randomization Approach](https://doi.org/10.1177/0962280219851817)" 
by Ioan Gabriel Bucur, Tom Claassen and Tom Heskes. The data set also contains 
simulated data necessary for exactly reproducing the figures in the article as
well as the routines necessary for recreating it. This research is presented 
in Chapter 4 of the PhD thesis titled "Being Bayesian about Causal Inference" by
Ioan Gabriel Bucur. 

The code is written in the R and C++ programming languages. BayesMR makes use of
the lite version of [PolyChord](https://github.com/PolyChord/PolyChordLite) 
nested sampling algorithm, which is owned and copyrighted by Will Handley, Mike 
Hobson and Anthony Lasenby. The PolyChordLite source code is bundled (as-is, 
except for build and install configuration) as a submodule in this package. 
For more details, please see the README and LICENSE accompanying the PolyChordLite 
submodule. This package also makes use of the [simpleini](https://github.com/brofield/simpleini) 
API for reading and writing INI configuration files, which is owned and copyrighted
by Brodie Thiesfield. For more details, please the README and LICENSE accompanying 
the simpleini submodule.


## Structure

The code is structured on the skeleton of an [R package](https://r-pkgs.org/index.html) 
package as follows:

- The folder `data` contains pre-saved simulated data in .RData format, which we 
use to recreate the figures from the article. The simulated data can also be
reproduced using the `reproduce-data.R` file from the `scripts` folder. 
The simulated data sets are described in `R/data.R`.

- The folder `inst` contains the package submodules.
  - The subfolder `inst/PolyChordLite` contains the [PolyChordLite](https://github.com/igbucur/PolyChordLite) submodule, which is
  a fork of the project developed by Will Handley, Mike Hobson and Anthony Lasenby.
  The PolyChordLite submodule implements a complex nested sampling algorithm,
  which we employ in BayesMR for estimating model evidences and producing posterior 
  samples. The original implementation is kept unchanged, with the exception of 
  changes in the build configuration file made for easier compilation and installation.
  - The subfolder `inst/simpleini` contains the [simpleini](https://github.com/igbucur/simpleini)
  submodule, which is a fork of the project developed by Brodie Thiesfield. 
  The simpleini submodule is a cross-platform library that provides a simple 
  API to read and write INI-style configuration files, which we use as command 
  line arguments for configuring each run of the BayesMR program.
  - The subfolder `inst/extdata` contains a few data files that are not in .RData
  format. `Parkinson_BMI_genetic_associations.csv` is a comma-separated values 
  file which contains a list of genetic variants and their associations with
  BMI and the risk of Parkinson's disease from [(Noyce et al., 2017)](https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1002314).
  `BayesMR_SS.txt` and `BayesMR_sigma_G.txt` are first-order and second-order statistics
  used as input for the provided BayesMR example (see also `ini/BayesMR.ini`).

- The folder `R` contains the R files necessary for reproducing the figures from
the article.

- The folder `man` contains the documentation for the implemented R functions
and for the registered .RData files.

- The folder `posterior` contains the C++ implementation of BayesMR written for
integration with PolyChord.

- The folder `ini` contains an example configuration file (`BayesMR.ini`). This
configuration file can be passed as the single necessary command line argument
when running BayesMR.

- The folder `scripts` contains the script `smmr-article-figures.R`, which
can be run from R to produce the figures in the Statistical Methods in Medical
Research (SMMR) article and a script called `reproduce-data.R`, which can be
used to reproduce the simulated data saved in the `data` folder.

- The top folder also contains the following files:
  - `DESCRIPTION` is the file describing the R package.
  - `NAMESPACE` is the file specifying the fucntions provided by the R package.
  - `LICENSE.md` is the file containing the GPL-3 license.


## Prerequisites

Just for installing the software and running basic examples, you will need:

- [GNU Bash](https://www.gnu.org/software/bash/)
- [GNU Make Build Tool](https://www.gnu.org/software/make/)
- [GNU C++ Compiler](https://gcc.gnu.org/) or similar
- [GNU Fortran Compiler](https://gcc.gnu.org/fortran/) or similar
- [Armadillo C++ Library for Linear Algebra](http://arma.sourceforge.net/)
- [Boost C++ Libraries](https://www.boost.org/)

These prerequisites are routinely available for Linux-based environments. For 
example, the (missing) prerequisites can be installed on Ubuntu using the following 
terminal command: `sudo apt install make lib{armadillo|boost}-dev`.

This software can also run on Windows, but typically requires quite a bit more 
work for installing all the prerequisites. The simplest way to make it work 
would be to install the most recent version of the [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10) (WSL). The
WSL can be used to emulate a Linux distribution like Ubuntu on a machine running 
Windows, so the prerequisites can be installed as if in a Linux-based environment.

A slightly more complicated but workable and native solution would be to employ 
the [MSYS2](https://www.msys2.org/) build platform. This setup involves installing
packages (including necessary dependencies) from the MSYS2 collection: `pacman -S mingw-w64-x86_64-{armadillo|boost|gcc-fortran|hdf5}`.

For installing and running the BayesMR R package, a few R packages are required. 
These are specified in the package `DESCRIPTION` file. These can also be
installed from R by running (from the BayesMR folder):
```
install.packages('devtools') # required package
devtools::install_deps(".", dependencies = TRUE) # install BayesMR package dependencies
```


## Installation Instructions

1. Download the software from GitHub with the following command:
`git clone --recurse-submodules https://github.com/igbucur/BayesMR.git`.

2. Enter the folder where the repository was cloned (called BayesMR by default).
Build the PolyChord nested sampling software and BayesMR by simply running `make` 
in the root directory. It is possible that some compiler flags will have to be set 
individually, so that the C++ compiler, e.g., [gcc](https://gcc.gnu.org/), can 
find the required headers and libraries. On Windows, the path to the collection
of prerequisites must be specified separately. For example, if using 
[MSYS2](https://www.msys2.org/) with 64-bit [MinGW](http://www.mingw.org/), then 
one must add `-I/c/msys64/mingw64/include` to `CPPFLAGS` and `-L/c/msys64/mingw64/lib` 
to `LDFLAGS` in the Makefile, assuming default installation directories.

Assuming everything managed to compile successfully, the executable `BayesMR` 
should appear in the root folder. To verify that BayesMR has been built successfully, 
one can use the provided script by running `bash scripts/run_BayesMR_example.sh` in a 
Linux or Windows terminal (requires [GNU Bash](https://www.gnu.org/software/bash/)). 
To verify that the R package is installed successfully, run `Rscript scripts/smmr-article-figures.R`.


## Licensing

BayesMR algorithm - A Bayesian Mendelian Randomization Approach

Copyright (C) 2020 Ioan Gabriel Bucur <ioan.gabriel.bucur@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
