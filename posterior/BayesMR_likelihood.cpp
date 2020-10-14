#include "BayesMR_likelihood.hpp"
#define _USE_MATH_DEFINES // necessary for M_PI on Windows
#include <cmath>

#include <iostream>
#include <boost/math/distributions/normal.hpp>
#include <armadillo>
using namespace arma;

#include "../inst/simpleini/SimpleIni.h" // For reading INI files

// This module is where your likelihood code should be placed.
//
// * The loglikelihood is called by the subroutine loglikelihood.
// * The likelihood is set up by setup_loglikelihood.
// * You can store any global/saved variables in the module.


//============================================================
// insert likelihood variables here
//
//
//============================================================

typedef struct {
  unsigned J; // number of instrumental variables
  mat SS;
  double slab_precision;
  double spike_precision;
  int model; // 0 - Slab, 1 - Spike, 2 - Reverse
  unsigned N; // number of observations
  vec sigma_G;
  boost::math::normal spike_gaussian;
  boost::math::normal slab_gaussian;
} MyConfig; 

MyConfig config;

// name of file containing model configuration
std::string model_config_filename;

/* GLOBALS */
boost::math::normal standard_gaussian(0.0, 1.0);
boost::math::normal uninformative_gaussian(0.0, 10.0);


double quantile_spike_and_slab(double F, double beta, 
    boost::math::normal& slab, boost::math::normal& spike) {
  
  if (fabs(beta - 1) < 1e-12) {
    return quantile(slab, F);
  } else if (fabs(beta) < 1e-12) {
    return quantile(spike, F);
  }
  
  double left = F < 0.5 ? quantile(slab, F) : quantile(spike, F);
  double left_value = beta * cdf(slab, left) + (1 - beta) * cdf(spike, left) - F;
  
  double right = F < 0.5 ? quantile(spike, F) : quantile(slab, F);
  double right_value = beta * cdf(slab, right) + (1 - beta) * cdf(spike, right) - F;
  
  double middle = (left + right) / 2;
  double mid_value = beta * cdf(slab, middle) + (1 - beta) * cdf(spike, middle) - F;
  
  while(fabs(mid_value) > 1e-12) {
    if (left_value * mid_value < 0) {
      right = middle;
      right_value = mid_value;
    } else {
      left = middle;
      left_value = mid_value;
    }
    middle = (left + right) / 2;
    mid_value = beta * cdf(slab, middle) + (1 - beta) * cdf(spike, middle) - F;
  }
  return middle;
}


// Main loglikelihood function
//
// Either write your likelihood code directly into this function, or call an
// external library from it. This should return the logarithm of the
// likelihood, i.e:
//
// loglikelihood = log_e ( P (data | parameters, model ) )
//
// theta are the values of the input parameters (in the physical space 
// NB: not the hypercube space).
//
// nDims is the size of the theta array, i.e. the dimensionality of the parameter space
//
// phi are any derived parameters that you would like to save with your
// likelihood.
//
// nDerived is the size of the phi array, i.e. the number of derived parameters
//
// The return value should be the loglikelihood
//
// This function is called from likelihoods/fortran_cpp_wrapper.f90
// If you would like to adjust the signature of this call, then you should adjust it there,
// as well as in likelihoods/my_cpp_likelihood.hpp
// 
double loglikelihood (double theta[], int nDims, double phi[], int nDerived)
{
  double logL = 0.0;

  vec salpha(config.J, fill::zeros);
  vec sgamma(config.J, fill::zeros);
  double sbeta = 0.0;
  double skappa_X = 0.0;
  double skappa_Y = 0.0;  

  unsigned par_no = 0L;
 
  // uniform prior on Gaussian mixture weights (indices 1:2)
  double wgamma = theta[par_no++];
  double walpha = theta[par_no++];
  
  // The default model assumes spike-and-slab priors for every structural parameter and pooled weights for the spike-and-slab components
  if (config.model == 0) {

    // Prior on pleiotropic effects (indices 2:(2J+1))
    for(int j = 0; j < config.J; ++j) {

      sgamma(j) = quantile_spike_and_slab(theta[par_no++], wgamma, config.slab_gaussian, config.spike_gaussian);
      salpha(j) = quantile_spike_and_slab(theta[par_no++], walpha, config.slab_gaussian, config.spike_gaussian);
    }

    // Prior on causal effect
    sbeta = quantile_spike_and_slab(theta[par_no++], 0.5, config.slab_gaussian, config.spike_gaussian);

    // Prior on confounding coefficients
    skappa_X = quantile_spike_and_slab(theta[par_no++] / 2 + 0.5, 0.5, config.slab_gaussian, config.spike_gaussian);
    skappa_Y = quantile_spike_and_slab(theta[par_no++], 0.5, config.slab_gaussian, config.spike_gaussian);

  // This model corresponds closely to the IV model assumptions (See Figure 18 in SMMR paper)
  } else if (config.model == 1) {

    // Prior on pleiotropic effects (indices 2:(2J+1))
    for(int j = 0; j < config.J; ++j) {
      sgamma(j) = quantile(config.slab_gaussian, theta[par_no++]); // Gaussian prior on instrument strengths
      salpha(j) = quantile(config.spike_gaussian, theta[par_no++]); // very small (irrelevant) pleiotropy
    }

    // Prior on causal effect
    sbeta = quantile(config.slab_gaussian, theta[par_no++]);

    // Prior on confounding coefficients
    skappa_X = quantile(config.slab_gaussian, theta[par_no++] / 2 + 0.5);
    skappa_Y = quantile(config.slab_gaussian, theta[par_no++]);
  }



  // Prior on scale parameters (half-normal)
  double sigma_X = quantile(uninformative_gaussian, theta[par_no++] / 2 + 0.5);
  double sigma_Y = quantile(uninformative_gaussian, theta[par_no++] / 2 + 0.5);

  assert(par_no == 2 * config.J + 7);
 
  mat f(config.J + 3, 2, fill::zeros); 
 
  // NOTE: not needed when mu_X = mu_Y = 0
  // f(0, 0) = - sigma_X * mu_X;
  // f(0, 1) = - sigma_Y * (mu_Y + sbeta * mu_X);

  for(int j = 0; j < config.J; ++j) {
#ifdef DEBUG
    std::cout << "(j, sigma_X, sgamma(j), config.sigma_G(j)): " << j << " " << sigma_X << " " << sgamma(j) << " " << config.sigma_G(j) << std::endl;
#endif
    f(j+1, 0) = - sigma_X * sgamma(j) / config.sigma_G(j);
    f(j+1, 1) = - sigma_Y * (salpha(j) + sbeta * sgamma(j)) / config.sigma_G(j);
  }
  f(config.J + 1, 0) = 1;
  f(config.J + 2, 1) = 1;
  
#ifdef DEBUG
  std::cout << "f: " << f << std::endl;
#endif

  mat S = f.t() * config.SS * f;

#ifdef DEBUG
  std::cout << "S: " << S << std::endl;
#endif

  mat Sigma(2, 2, fill::zeros);
  Sigma(0, 0) = sigma_X * sigma_X * (1 + skappa_X * skappa_X);
  Sigma(0, 1) = Sigma(1, 0) = sigma_X * sigma_Y * (skappa_X * skappa_Y + sbeta * (1 + skappa_X * skappa_X));
  Sigma(1, 1) = sigma_Y * sigma_Y * (1 + sbeta * sbeta + (skappa_Y + sbeta * skappa_X) * (skappa_Y + sbeta * skappa_X));

#ifdef DEBUG
  std::cout << "sigma: " << Sigma << std::endl;
#endif

  logL = config.N * (-0.5 * (log(4 * M_PI * M_PI * det(Sigma)) + trace(Sigma.i() * S)));

#ifdef DEBUG
  std::cout << "end call ll: " << logL << std::endl;
#endif

  // NOTE: not needed when theta is fixed
  // + mean_G.t() * log(theta) + (config.n - mean_G).t() * log(1 - theta);  
    
  return logL;

}

// Prior function
//
// Either write your prior code directly into this function, or call an
// external library from it. This should transform a coordinate in the unit hypercube
// stored in cube (of size nDims) to a coordinate in the physical system stored in theta
//
// This function is called from likelihoods/fortran_cpp_wrapper.f90
// If you would like to adjust the signature of this call, then you should adjust it there,
// as well as in likelihoods/my_cpp_likelihood.hpp
// 
void prior (double cube[], double theta[], int nDims)
{
  
  for (int i = 0; i < nDims; ++i) {
    theta[i] = cube[i];
  }  
}

// Dumper function
//
// This function gives you runtime access to variables, every time the live
// points are compressed by a factor settings.compression_factor.
//
// To use the arrays, subscript by following this example:
//
//    for (auto i_dead=0;i_dead<ndead;i_dead++)
//    {
//        for (auto j_par=0;j_par<npars;j_par++)
//            std::cout << dead[npars*i_dead+j_par] << " ";
//        std::cout << std::endl;
//    }
//
// in the live and dead arrays, the rows contain the physical and derived
// parameters for each point, followed by the birth contour, then the
// loglikelihood contour
//
// logweights are posterior weights
// 
void dumper(int ndead,int nlive,int npars,double* live,double* dead,double* logweights,double logZ, double logZerr)
{
}


// Ini path reading function
//
// If you want the file path to the ini file (for example for pointing to other config files), store ite value in a global variable with this function
void set_ini(std::string ini_str_in)
{
    /**
    Set value for constant holding ini file path.

    @param ini_str_in: ini file path
    **/
  model_config_filename = ini_str_in;
}


// Setup of the loglikelihood
// 
// This is called before nested sampling, but after the priors and settings
// have been set up.
// 
// This is the time at which you should load any files that the likelihoods
// need, and do any initial calculations.
// 
// This module can be used to save variables in between calls
// (at the top of the file).
// 
// All MPI threads will call this function simultaneously, but you may need
// to use mpi utilities to synchronise them. This should be done through the
// integer mpi_communicator (which is normally MPI_COMM_WORLD).
//
// This function is called from likelihoods/fortran_cpp_wrapper.f90
// If you would like to adjust the signature of this call, then you should adjust it there,
// as well as in likelihoods/my_cpp_likelihood.hpp
//
void setup_loglikelihood()
{

  CSimpleIniA ini;
  ini.LoadFile(model_config_filename.c_str());

  std::string SS_filename = ini.GetValue("model settings", "SS_filename", "");
  std::string sigma_G_filename = ini.GetValue("model settings", "sigma_G_filename", "");
  config.slab_precision = strtod(ini.GetValue("model settings", "slab_precision", ""), NULL); 
  config.spike_precision = strtod(ini.GetValue("model settings", "spike_precision", ""), NULL); 
  config.J = strtoul(ini.GetValue("model settings", "instruments", ""), NULL, 0);
  config.N = strtoul(ini.GetValue("model settings", "observations", ""), NULL, 0);
  config.model = strtol(ini.GetValue("model settings", "model", ""), NULL, 0);
 
  std::cout << "SS filename: " << SS_filename << std::endl;
  
  config.SS.load(SS_filename);
  std::cout << "SS matrix: " << std::endl;
  std::cout << config.SS << std::endl;

  std::cout << "sigma_G filename: " << sigma_G_filename << std::endl;
  
  config.sigma_G.load(sigma_G_filename);
  std::cout << "sigma_G matrix: " << std::endl;
  std::cout << config.sigma_G << std::endl;

  std::cout << "Slab precision: " << config.slab_precision << std::endl;
  std::cout << "Spike precision: " << config.spike_precision << std::endl;
  std::cout << "Number of instruments: " << config.J << std::endl;
  std::cout << "Number of samples: " << config.N << std::endl;
  std::cout << "Model: " << config.model << std::endl;


  config.spike_gaussian = boost::math::normal(0.0, sqrt(1.0 / config.spike_precision));
  config.slab_gaussian = boost::math::normal(0.0, sqrt(1.0 / config.slab_precision));

}
