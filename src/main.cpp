/**
 * Copyright (C) 2021 Ante Qu <antequ@cs.stanford.edu>.
 */

#include "define_args.h"
#include "linking_number_certificate.h"
#include "model.h"

#include <chrono>
#include <iostream>

using namespace VerifyCurves;

int main(int argc, char **argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  Model model;
  // Import the input model
  {
    std::cout << "Importing from " << FLAGS_input << "." << std::endl;
    std::chrono::steady_clock::time_point import_begin =
        std::chrono::steady_clock::now();
    model.ImportFromFile(FLAGS_input);

    std::chrono::steady_clock::time_point import_end =
        std::chrono::steady_clock::now();
    int64_t import_time = std::chrono::duration_cast<std::chrono::microseconds>(
                              import_end - import_begin)
                              .count();
    std::cout << "Import runtime: " << (1e-3 * import_time) << " ms."
              << std::endl;
  }

  bool use_direct_sum = false;
  if (FLAGS_method.compare("directsum") == 0) {
    use_direct_sum = true;
  } else if (FLAGS_method.compare("barneshut") == 0) {
    use_direct_sum = false;
  } else {
    std::cerr << "Invalid method: '" << FLAGS_method << "'." << std::endl;
    throw std::runtime_error("Invalid method encountered.");
  }
  LinkingNumberCertificate certificate;

  // Time and compute the linking-number certificate.
  {
    std::string method_string =
        use_direct_sum ? "Direct Summation" : "Barnes–Hut";
    std::cout << "Computing the certificate using the " << method_string
              << " Method." << std::endl;
    std::chrono::steady_clock::time_point compute_begin =
        std::chrono::steady_clock::now();
    certificate.ComputeFromModel(model, use_direct_sum, FLAGS_barnes_hut_init_beta,
                                 FLAGS_barnes_hut_beta_limit);
    std::chrono::steady_clock::time_point compute_end =
        std::chrono::steady_clock::now();
    int64_t compute_time =
        std::chrono::duration_cast<std::chrono::microseconds>(compute_end -
                                                              compute_begin)
            .count();
    std::cout << "Computation runtime: " << (1e-3 * compute_time)
              << " ms. This includes all 3 stages." << std::endl;
  }

  // Time and export the linking-number certificate.
  {
    std::cout << "Exporting the certificate to " << FLAGS_output << "."
              << std::endl;
    std::chrono::steady_clock::time_point export_begin =
        std::chrono::steady_clock::now();
    certificate.ExportToFile(FLAGS_output);
    std::chrono::steady_clock::time_point export_end =
        std::chrono::steady_clock::now();
    int64_t export_time = std::chrono::duration_cast<std::chrono::microseconds>(
                              export_end - export_begin)
                              .count();
    std::cout << "Export runtime: " << (1e-3 * export_time) << " ms."
              << std::endl;
  }

  // Perform the comparison if the user requests it.
  if (FLAGS_comparewith.compare("[none]") != 0) {
    if (FLAGS_diffpairs.compare("[none]") != 0 ||
        FLAGS_diffcurves.compare("[none]") != 0) {
      LinkingNumberCertificate reference_certificate(FLAGS_comparewith);
      certificate.ExportDiffsToFiles(FLAGS_diffpairs, FLAGS_diffcurves,
                                     reference_certificate);
      std::cout
          << "Compared and exported the difference between the model and the "
             "reference certificate."
          << std::endl;
    }
  }
  return 0;
}
