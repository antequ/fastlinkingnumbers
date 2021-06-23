/**
 * Copyright (C) 2021 Ante Qu <antequ@cs.stanford.edu>.
 */

#include "model.h"
#include <chrono>
#include <gflags/gflags.h>
#include <iostream>

DEFINE_string(
    input, "data/houdiniobj.obj",
    "The filepath to the OBJ input file, which must end with '.obj'.");

DEFINE_string(
    output, "data/houdiniobj.bcc",
    "The filepath to the BCC output file, which must end with '.bcc'.");

DEFINE_string(
    segment_type, "PL",
    "The desired segment type for the BCC output file, which must be either "
    "uniform Catmull-Rom (C0), B-Spline (BS), Bezier (BZ), or polyline (PL). "
    "You must specify it with the two-character abbreviation, which is "
    "case-sensitive.");

using namespace VerifyCurves;

int main(int argc, char **argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  Model model;
  SegmentType segment_type;
  if (FLAGS_segment_type.compare("C0") == 0) {
    segment_type = SegmentType::CatmullRom;
  } else if (FLAGS_segment_type.compare("BS") == 0) {
    segment_type = SegmentType::BSpline;
  } else if (FLAGS_segment_type.compare("BZ") == 0) {
    segment_type = SegmentType::Bezier;
  } else if (FLAGS_segment_type.compare("PL") == 0) {
    segment_type = SegmentType::Polyline;
  } else {
    throw std::runtime_error(
        "Invalid segment_type. (Must be 'C0', 'BS', 'BZ', or 'PL'.)");
  }

  // OBJ Import
  {
    std::cout << "Importing from " << FLAGS_input << "." << std::endl;
    std::chrono::steady_clock::time_point import_begin =
        std::chrono::steady_clock::now();
    model.ImportFromObjFile(FLAGS_input, segment_type);

    std::chrono::steady_clock::time_point import_end =
        std::chrono::steady_clock::now();
    int64_t import_time = std::chrono::duration_cast<std::chrono::microseconds>(
                              import_end - import_begin)
                              .count();
    std::cout << "Import runtime: " << (1e-3 * import_time) << " ms."
              << std::endl;
  }

  // BCC Export
  {
    std::cout << "Exporting to " << FLAGS_output << "." << std::endl;
    std::chrono::steady_clock::time_point export_begin =
        std::chrono::steady_clock::now();
    model.ExportToFile(FLAGS_output);
    std::chrono::steady_clock::time_point export_end =
        std::chrono::steady_clock::now();
    int64_t export_time = std::chrono::duration_cast<std::chrono::microseconds>(
                              export_end - export_begin)
                              .count();
    std::cout << "Export runtime: " << (1e-3 * export_time) << " ms."
              << std::endl;
  }

  return 0;
}
