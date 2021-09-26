/**
 * Copyright (C) 2021 Ante Qu <antequ@cs.stanford.edu>.
 */

#include "define_args.h"
#include "linking_number_certificate.h"
#include "model.h"

#include <fstream>
#include <iostream>

using namespace VerifyCurves;

// A table1 row consists of:
// Model Name, N_I, L, P, PLS time, Discr. time, DS time, BH time, BH (average)
// error A table2 row consists of: Model Name, N_I, N_D
void PrintTableRow(std::ostream &out1, std::ostream &out2,
                   const std::string &filename,
                   const std::string &referencefilename,
                   const std::string &modelname) {
  Model model;
  // Import the input model
  model.ImportFromFile(filename);
  // Compute N_I, L
  int N_I = 0;
  int L = model.GetCurves().size();
  for (int i = 0; i < L; ++i) {
    N_I += model.GetCurves()[i].num_segments();
  }
  // Load reference matrix
  LinkingNumberCertificate certificate;
  Eigen::SparseMatrix<int> reference =
      certificate.ImportFromFile(referencefilename);
  // Time and compute the linking-number certificate.
  double avg_computed_error;
  double max_computed_error;
  int N_D, P;
  double pls_time, discr_time, ds_time, bh_time;
  bool compute_both = true;
  // disable compute_both for rows 13, 14, 16, 17, 18, 19, 20, because Direct
  // Sum takes too long.
  if ((filename.find("N20M") != std::string::npos) ||
      (filename.find("N4M") != std::string::npos) ||
      (filename.find("N40M") != std::string::npos) ||
      (filename.find("N1M") != std::string::npos) ||
      (filename.find("N2M") != std::string::npos)) {
    std::cout << "Disabling direct sum compute for " << modelname << std::endl;
    compute_both = false;
  }

  ds_time = -1;
  certificate.ComputeFromModelReplicability(
      model, true /* compare_against_reference */, reference,
      avg_computed_error, max_computed_error, N_D, P, pls_time, discr_time,
      ds_time, bh_time, compute_both /* compute_both */,
      false /* force_direct_sum */, FLAGS_barnes_hut_init_beta,
      FLAGS_barnes_hut_beta_limit);
  out1 << "\"" << modelname << "\", " << N_I << ", " << L << ", " << P << ", "
       << pls_time << ", " << discr_time << ", "
       << ((ds_time < -0.5) ? "N/A" : std::to_string(ds_time)) << ", "
       << bh_time << ", " << avg_computed_error << std::endl;
  std::cout << "Model " << modelname << ", Avg error: " << avg_computed_error
            << ", Max error: " << max_computed_error << std::endl;
  out2 << "\"" << modelname << "\", " << N_I << ", " << N_D << std::endl;
}

int main(int argc, char **argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  std::cout << "This script will generate Tables 1 and 2 from the paper, into "
               "table1.csv and table2.csv, respectively."
            << std::endl;
  std::cout << "Please ensure that the dataset is in the ./dataset/ folder."
            << std::endl;
  std::ofstream table1("table1.csv", std::ofstream::out);
  std::ofstream table2("table2.csv", std::ofstream::out);
  std::string dataset_path = "dataset/";
  std::string reference_path = dataset_path + "reference_certificates/";
  // Print table headers
  // A table1 row consists of:
  // Model Name, N_I, L, P, PLS time, Discr. time, DS time, BH time, BH
  // (average) error
  table1 << "Model, N_I, L, P, PLS Time, Dscr. Time, DS Compute Time, BH "
            "Compute Time, Avg. BH Abs. Error"
         << std::endl;
  // A table2 row consists of:
  // Model Name, N_I, N_D
  table2 << "Model, N_I, N_D" << std::endl;
  std::vector<std::string> filenames = {"alien_sweater_init",
                                        "alien_sweater_final",
                                        "sheep_sweater",
                                        "sweater",
                                        "glove",
                                        "knittubeinit",
                                        "knittubefinal",
                                        "chainmail_init",
                                        "chainmail_final",
                                        "rubber_bands_final",
                                        "double_helix_ribbon_lambda10_N200K",
                                        "double_helix_ribbon_lambda1K_N200K",
                                        "double_helix_ribbon_lambda10_N20M",
                                        "double_helix_ribbon_lambda1K_N20M",
                                        "thicksquarelink_N500K",
                                        "thicksquarelink_N4M",
                                        "torus_lambda1M_N20M",
                                        "torus_lambda0_N40M",
                                        "woundball_nu1K_N1M",
                                        "woundball_nu10K_N2M"};
  std::vector<std::string> modelnames = {"Alien Sweater (Initial)",
                                         "Alien Sweater (Final)",
                                         "Sheep Sweater",
                                         "Sweater",
                                         "Glove",
                                         "Knit Tube (Initial)",
                                         "Knit Tube (Final)",
                                         "Chainmail (Initial)",
                                         "Chainmail (Final)",
                                         "Rubber Bands",
                                         "Double-Helix Ribbon λ=10, 200K Segs",
                                         "Double-Helix Ribbon λ=1K, 200K Segs",
                                         "Double-Helix Ribbon λ=10, 20M Segs",
                                         "Double-Helix Ribbon λ=1K, 20M Segs",
                                         "Thick Square Link, 500K Segs",
                                         "Thick Square Link, 4M Segs",
                                         "Torus λ=1M, 20M Segs",
                                         "Torus λ=0,  40M Segs",
                                         "Woundball ν=1K,  1M Segs",
                                         "Woundball ν=10K, 2M Segs"};
  assert(filenames.size() == modelnames.size());
  for (int i = 0; i < filenames.size(); ++i) {
    std::string filename = dataset_path + filenames[i] + ".bcc";
    std::string referencefilename = reference_path + filenames[i] + ".txt";
    const std::string &modelname = modelnames[i];
    PrintTableRow(table1, table2, filename, referencefilename, modelname);
  }
  table1.close();
  table2.close();
  return 0;
}
