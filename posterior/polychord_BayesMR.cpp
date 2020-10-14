#include "../inst/PolyChordLite/src/polychord/interfaces.hpp"
#include "BayesMR_likelihood.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    if (argc == 2) {
        std::string config_file = argv[1];

        set_ini(config_file); // pass config file to likelihood computation
        run_polychord(loglikelihood,setup_loglikelihood, config_file) ;
        return 0;
    }
    else{
        std::cerr << "BayesMR must be called with exactly one file specifying a INI configuration file. The format of the configuration file is given in ini/BayesMR.ini. Configuration files can be easily created using the R function R/write_BayesMR_configuration_file.R provided in the BayesMR package." << std::endl;
        return 1;
    }
}

