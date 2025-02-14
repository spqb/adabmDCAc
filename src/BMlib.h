// BOLTZMANN MACHINE CODE - version of March, 2021
// Authors: Anna Paola Muntoni and Francesco Zamponi
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <getopt.h>
#include <stdbool.h>
#include <valarray>
#include <iostream>
#include "BMaux.h"

using namespace std;
typedef double MYFLOAT;

#ifndef BMLIB_H
#define BMLIB_H

class Params
{
public:
  FILE * filel;
  char *file_msa, *file_freq, *file_w, *file_params, *file_msa_e, init, *label, *outputfolder, *ctype, *file_3points, *file_cc, *file_last_chain;
  bool dms_flag, energies_flag, restore_flag, deczero, print_samples, Metropolis, Gibbs, nprinteq, rmgauge, dgap, gapnn, phmm, compwise, persistent, initdata, overwrite, adapt, dec_sdkl, get_abc;
  double initst, sparsity, rho, w_th, regJ1, regJ2, lrateJ, lrateh, conv, pseudocount, betaJ, betaH, lambda_e, drate, nactive;
  int tau, seed, learn_strat, nprint, nprintfile, Teq, Nmc_starts, Nmc_config, Twait, Twait_last, maxiter, gsteps, nmix, num_threads;

  Params();

  int read_params(int &argc, char **argv);
  void print_learning_strategy();
  void construct_filenames(int iter, bool conv, char *par, char *par_zsum, char *ene, char *corr, char *score, char *first, char *sec, char *third, char *lchain, char *eqfile);
};

class Data
{
public:
  int q, L, M;
  double Meff;
  vector<vector<unsigned char> > msa;
  vector<MYFLOAT> w;
  vector<MYFLOAT> fm;
  vector<vector<MYFLOAT> > sm;
  vector<vector<MYFLOAT> > cov;
  vector<MYFLOAT> tm; // Pay attention. File contains connected 3rd order correlations
  vector<vector<int> > tm_index;
  Params *params;

  Data(Params *_params);

  void read_alphabet();
  void read_msa();
  void compute_w();
  void alloc_structures();
  void compute_empirical_statistics();
  void read_freq();

  /******************** METHODS FOR 3RD ORDER STATISTICS ***********************************************************/

  void load_third_order_indices();
  /******************** METHODS FOR OUTPUT ***********************************************************/

  void print_msa(char *filename);
  int print_statistics(char *file_sm, char *file_fm, char *file_tm, char *file_c, vector<MYFLOAT> &fm_s, vector<vector<MYFLOAT> > &sm_s, vector<MYFLOAT> &tm_s);
};

/* CONTAINER FOR THE ENERGY MSA AND ITS FITNESS */
class Data_e
{
public:
  Params *params;
  int q, L, M;
  double A;
  vector<vector<unsigned char> > msa;
  vector<double> w;
  vector<MYFLOAT> fitness;
  vector<MYFLOAT> energy;
  vector<MYFLOAT> tg_energy;
  vector<MYFLOAT> gradh;
  vector<vector<MYFLOAT> > gradJ;

  Data_e(Params *_params, int _q, int _L, double Meff);
  void read_msa();
  void compute_w();
  void compute_grad();
  void set_tg_energy();
  void set_tg_energy_rank();
  void set_tg_energy_fitness();
  double demeanwt();
  double demean();
  double destd();

  void print_msa(char *filename);
  void print_energy(char *filename);
};
#endif
