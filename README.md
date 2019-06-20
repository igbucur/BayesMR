# BayesMR

Code implementation for the approach described in the article "Inferring the Direction of a Causal Link and Estimating Its Effect via a Bayesian Mendelian Randomization Approach" by Ioan Gabriel Bucur, Tom Claassen, and Tom Heskes (https://doi.org/10.1177/0962280219851817).

## PolyChord
The BayesMR approach uses the PolyChord software (Copyright © 2015 Will Handley, Mike Hobson & Anthony Lasenby) forked from the GitHub repository https://github.com/PolyChord/PolyChordLite. The license agreement can be found at https://github.com/PolyChord/PolyChordLite/blob/master/LICENCE or within the PolyChordLite submodule included in this repository.

## Installation instructions
1. Clone this repository using the --recurse-submodules option, e.g.:
`git clone --recurse-submodules https://github.com/igbucur/BayesMR.git`.

2. Enter the folder where the repository was cloned (called BayesMR by default) and execute the Makefile using: `make polychord_MR`. Assuming everything managed to compile successfully, the PolyChord implementation of the BayesMR approach should appear as the executable file _BayesMR/bin/polychord_MR_.

3. Execute the script `./run_SMMR_experiments.sh` provided with the package to verify that the installation is working.
