/**
 * Copyright (C) 2021 Ante Qu <antequ@cs.stanford.edu>.
 * Copyright (C) 2018 Jonathan Leaf.  All rights reserved.
 * Originally "sim_model.cpp" by Jonathan Leaf. Modified by Ante Qu.
 */

#include "model.h"

#include "curve.h" // for Point3, SegmentType, and Curve.

#include <algorithm> // for transform
#include <cctype>    // for tolower
#include <fstream>
#include <iostream>
#include <string>

namespace VerifyCurves {
/**
 * The default constructor generates an empty list of curves.
 */
Model::Model() {}

// This helper function converts strings to lowercase, which we use in checking
// file extensions. It's from
// https://en.cppreference.com/w/cpp/string/byte/tolower.
std::string str_tolower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 // static_cast<int(*)(int)>(std::tolower)         // wrong
                 // [](int c){ return std::tolower(c); }           // wrong
                 // [](char c){ return std::tolower(c); }          // wrong
                 [](unsigned char c) { return std::tolower(c); } // correct
  );
  return s;
}

/**
 * This method imports a curve collection from a BCC file and stores it in
 * this Model, in curves_.
 *
 * @param [in] filename the filepath of the BCC file to import from.
 *
 * @post curves_ is updated with the contents from the BCC file.
 */
void Model::ImportFromFile(const std::string &filename) {
  std::string fileext = str_tolower(filename.substr(filename.rfind(".")));
  int bcc_ext = strcmp(fileext.c_str(), ".bcc");
  if (bcc_ext == 0) {
    ImportBccFile(filename);
    segment_type_ = GetCurves()[0].get_segment_type();
  } else {
    std::cerr << "File extension for import must be .bcc: " << filename
              << std::endl;
    throw std::runtime_error("Input model file extension error.");
  }
}

/**
 * This method exports the curve collection in this Model to a BCC file.
 *
 * @param [in] filename the filepath of the BCC file to save to.
 */
void Model::ExportToFile(const std::string &filename) const {
  std::string fileext = str_tolower(filename.substr(filename.rfind(".")));
  int bcc_ext = strcmp(fileext.c_str(), ".bcc");
  if (bcc_ext == 0) {
    ExportBccFile(filename);
  } else {
    std::cerr << "File extension for export must be .bcc: " << filename
              << std::endl;
    throw std::runtime_error("Output model file extension error.");
  }
}

// This helper function returns the last character of ss. It's used in parsing
// OBJ files.
char peek_last_char(std::stringstream &ss) {
  // based on https://stackoverflow.com/questions/36743452/
  int last_char;
  std::stringstream::pos_type pos = ss.tellg();
  ss.seekg(-1, std::ios::end);
  last_char = ss.peek();
  ss.seekg(pos);
  if (last_char == EOF) {
    throw std::runtime_error("peek_last_char(): Line is empty.");
  }
  return last_char;
}

// This helper function removes leading and trailing whitespace from str. It's
// used in parsing OBJ files.
std::string trimfnc(const std::string &str) {
  // based on https://www.thecrazyprogrammer.com/2021/01/c-string-trim.html
  const char *typeOfWhitespaces = " \t\v";
  std::string::size_type endposition =
      str.find_last_not_of(typeOfWhitespaces) + 1;
  std::string::size_type startposition =
      str.find_first_not_of(typeOfWhitespaces);
  if (endposition < startposition)
    return "";
  else
    return str.substr(startposition, endposition - startposition);
}

// This helper function advances ss until the next character is not a
// whitespace. It's used in parsing OBJ files.
void eatspaces(std::stringstream &ss) {
  int next = ss.peek();
  while ((next != EOF) && (next == ' ' || next == '\t' || next == '\v')) {
    ss.get();
    next = ss.peek();
  }
}

/**
 * This method imports a curve collection from an OBJ file and stores it in
 * this Model, in curves_. It reads vertices from "v" and control points from
 * "l". It will attach adjacent lines ("l") together into the same curve if
 * their last and first vertices respectively match, which is necessary for
 * many Blender OBJs. All curves will be marked as closed loops.
 *
 * @param [in] filename the filepath of the OBJ file to import from.
 *
 * @post curves_ is updated with the contents from the BCC file.
 */
void Model::ImportFromObjFile(const std::string &filename,
                              SegmentType segment_type) {
  std::string fileext = str_tolower(filename.substr(filename.rfind(".")));
  int obj_ext = strcmp(fileext.c_str(), ".obj");
  if (obj_ext != 0) {
    std::cerr << "File extension must be OBJ: " << filename << std::endl;
    throw std::runtime_error("Import OBJ model: file extension error.");
  }
  // set up field types
  std::vector<Point3> vertices;
  int vid = -1;
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Unable to open file: " << filename << std::endl;
    throw std::runtime_error("Import OBJ model: unable to open file.");
  }
  std::string line;

  int linecount = 0;
  Curve curve;
  std::vector<Point3> &q = curve.get_q();
  int first_vid = -1;
  int last_vid = -1;
  curve.set_segment_type(segment_type);
  segment_type_ = segment_type;
  while (true) {
    if (!std::getline(file, line))
      break;
    ++linecount;
    if (linecount % 100000 == 0) {
      std::cout << "Processed the first " << linecount << " lines of "
                << filename << "." << std::endl;
    }
    line = trimfnc(line);
    if (line.length() == 0)
      continue;
    // ignore comment lines
    if (line.front() == '#')
      continue;
    // Houdini likes to output a long line for each curve. Attach lines split
    // apart by '\\'
    std::stringstream linestream;
    linestream << line;
    while (peek_last_char(linestream) == '\\') {
      if (!std::getline(file, line))
        break;
      line = trimfnc(line);
      std::stringstream::pos_type curr_pos = linestream.tellp();
      linestream.seekp(curr_pos - static_cast<std::streamoff>(1));
      // Note that this will append an extra space since the character before
      // '\\' is often a space.
      linestream << ' ' << line;
    }
    char first_char;
    linestream.get(first_char);
    int second_char = linestream.peek();
    if (second_char == EOF)
      continue;
    // We only parse "v " and "l " lines, so we require the second char to be a
    // whitespace.
    if (!((second_char == ' ') || (second_char == '\t') ||
          (second_char == '\v')))
      continue;
    if (first_char == 'v') {
      std::string value;
      // vertices
      double vx, vy, vz;
      bool successful = false;
      eatspaces(linestream);
      if (std::getline(linestream, value, ' ')) {
        vx = std::stod(value);
        eatspaces(linestream);
        if (std::getline(linestream, value, ' ')) {
          vy = std::stod(value);
          eatspaces(linestream);
          if (std::getline(linestream, value, ' ')) {
            vz = std::stod(value);
            successful = true;
          }
        }
      }
      if (!successful) {
        std::cerr << "Obj input: Invalid vertex (vid " << vertices.size() + 1
                  << ")." << std::endl;
        throw std::runtime_error("Obj input: Invalid vertex.");
      }
      vertices.push_back(Point3(vx, vy, vz));
    } else if (first_char == 'l') {
      std::string value;
      // polyline
      eatspaces(linestream);
      vid = -1;
      // Blender splits its curves onto multiple lines as line segments, so we
      // have to attach curves between lines and only finalize a curve once a
      // vertex is not shared.
      while (std::getline(linestream, value, ' ')) {
        bool first = (vid < 0);
        vid = std::stoi(value.substr(0, value.find('/')));
        if (first) {
          if (vid != last_vid) {
            // new curve
            if (q.size() > 1) {
              // Houdini repeats the first and last vertex, but Blender does
              // not.
              if ((last_vid == first_vid) || q.back().isApprox(q.front())) {
                q.pop_back();
              }
              curve.set_looped(true);
              curves_.push_back(curve);
              q.clear();
              last_vid = -1;
            }
            q.push_back(vertices[vid - 1]);
            last_vid = vid;
            first_vid = vid;
          }
        } else {
          if ((vid >= 1) && (vid < vertices.size() + 1)) {
            q.push_back(vertices[vid - 1]);
            last_vid = vid;
          } else {
            std::cerr << "Obj import: Illegal vid: " << vid << std::endl;
            throw std::runtime_error("Obj Import: illegal vertex ID.");
          }
        }
        eatspaces(linestream);
      }
    }
  }
  // insert the last curve
  if (q.size() > 1) {
    // Houdini repeats the first and last vertex, but Blender does not.
    if ((last_vid == first_vid) || q.back().isApprox(q.front())) {
      q.pop_back();
    }
    curve.set_looped(true);
    curves_.push_back(curve);
  }
  file.close();
  std::cout << "Successfully imported " << curves_.size() << " curves."
            << std::endl;
}

// This is a struct for loading BCC file headers, supplied from:
// http://www.cemyuksel.com/cyCodeBase/soln/using_bcc_files.html.
struct BCCHeader {
  char sign[3];
  unsigned char byteCount;
  char curveType[2];
  char dimensions;
  char upDimension;
  uint64_t curveCount;
  uint64_t totalControlPointCount;
  char fileInfo[40];
};

/// This method imports from a BCC file.
void Model::ImportBccFile(const std::string &filename) {
  // This code for loading BCC files is based on
  // http://www.cemyuksel.com/cyCodeBase/soln/using_bcc_files.html.
  BCCHeader header;
  FILE *pFile = fopen(filename.c_str(), "rb");
  if (pFile == nullptr) {
    std::cerr << "Error opening input file: " << filename << "\n" << std::endl;
    throw std::runtime_error("Unable to load BCC file.");
  }
  size_t readsize = fread(&header, sizeof(header), 1, pFile);
  if (readsize != 1) {
    std::cerr << "Incomplete header: " << filename << "\n" << std::endl;
    throw std::runtime_error("Unable to load BCC file.");
  }
  if ((header.sign[0] != 'B') || (header.sign[1] != 'C') ||
      (header.sign[2] != 'C')) {
    std::cerr << "Invalid header: " << filename << "\n" << std::endl;
    throw std::runtime_error("Unable to load BCC file.");
  }
  if (header.byteCount != 0x44) {
    std::cerr << "Invalid bytecount: " << filename
              << ". We only support 32-bit integers and floats.\n"
              << std::endl;
    throw std::runtime_error("Unable to load BCC file.");
  }
  if (header.dimensions != 3) {
    std::cerr << "Invalid dimensions: " << filename
              << ". We only support data in 3D Cartesian space.\n"
              << std::endl;
    throw std::runtime_error("Unable to load BCC file.");
  }
  // Read the curve type. It is:
  // C0 for uniform Catmull-Rom curves, BS for B-Splines,
  // BZ for Bezier curves, and PL for polylines.
  SegmentType curve_type;
  if ((header.curveType[0] == 'C') && (header.curveType[1] == '0')) {
    curve_type = SegmentType::CatmullRom;
  } else if ((header.curveType[0] == 'B') && (header.curveType[1] == 'S')) {
    curve_type = SegmentType::BSpline;
  } else if ((header.curveType[0] == 'B') && (header.curveType[1] == 'Z')) {
    curve_type = SegmentType::Bezier;
  } else if ((header.curveType[0] == 'P') && (header.curveType[1] == 'L')) {
    curve_type = SegmentType::Polyline;
  } else {
    std::cerr << "Invalid valid curve segment type: " << filename
              << ". We only support C0, BS, BZ, and PL.\n"
              << std::endl;
    throw std::runtime_error("Unable to load BCC file.");
  }
  // float buffer to be converted to doubles
  Eigen::Matrix3Xf point_buffer;
  curves_ = std::vector<Curve>(header.curveCount);
  for (int i = 0; i < curves_.size(); ++i) {
    int32_t control_point_count;
    readsize = fread(&control_point_count, sizeof(int32_t), 1, pFile);
    if (readsize != 1) {
      std::cerr << "Cannot read curve control point count: " << filename
                << ", curve number " << i << ".\n"
                << std::endl;
      throw std::runtime_error("Unable to load BCC file.");
    }
    curves_[i].set_looped(control_point_count < 0);
    curves_[i].set_segment_type(curve_type);
    if (control_point_count < 0) {
      control_point_count = -control_point_count;
    }
    if (control_point_count != 0) {
      point_buffer.resize(3, control_point_count);
      readsize = fread(point_buffer.data(), sizeof(float),
                       3 * control_point_count, pFile);
      if (readsize != 3 * control_point_count) {
        std::cerr << "Invalid number of control points read: " << filename
                  << ", curve number " << i << ", read " << readsize << " of "
                  << control_point_count << ".\n"
                  << std::endl;
        throw std::runtime_error("Unable to load BCC file.");
      }
      curves_[i].get_q().resize(control_point_count);
      Eigen::Map<Eigen::Matrix3Xd> casted_points(curves_[i].get_q()[0].data(),
                                                 3, control_point_count);
      casted_points = point_buffer.cast<double>();
    }
  }
  fclose(pFile);
}

/// This method exports to a BCC file.
void Model::ExportBccFile(const std::string &filename) const {
  // This code for writing out BCC files is based on
  // http://www.cemyuksel.com/cyCodeBase/soln/using_bcc_files.html.
  BCCHeader header;
  FILE *fp = fopen(filename.c_str(), "wb");
  if (fp == nullptr) {
    std::cerr << "Error creating output file: " << filename << std::endl;
    throw std::runtime_error("Unable to export BCC file.");
  }
  // Set header
  header.sign[0] = 'B';
  header.sign[1] = 'C';
  header.sign[2] = 'C';

  header.byteCount = 0x44; // 4 bytes for int, 4 bytes for float
  // Set the curve type. It is:
  // C0 for uniform Catmull-Rom curves, BS for B-Splines,
  // BZ for Bezier curves, and PL for polylines.
  // The default is C0 if there are no curves.
  header.curveType[0] = 'C';
  header.curveType[1] = '0';
  if (curves_.size() > 0) {
    if (curves_[0].get_segment_type() == SegmentType::CatmullRom) {
      header.curveType[0] = 'C';
      header.curveType[1] = '0';
    } else if (curves_[0].get_segment_type() == SegmentType::BSpline) {
      header.curveType[0] = 'B';
      header.curveType[1] = 'S';
    } else if (curves_[0].get_segment_type() == SegmentType::Bezier) {
      header.curveType[0] = 'B';
      header.curveType[1] = 'Z';
    } else if (curves_[0].get_segment_type() == SegmentType::Polyline) {
      header.curveType[0] = 'P';
      header.curveType[1] = 'L';
    }
  }
  header.dimensions = 3;
  // Set Y as the up dimension.
  header.upDimension = 1;
  header.curveCount = curves_.size();
  // compute total control point count
  header.totalControlPointCount = 0;
  for (int i = 0; i < header.curveCount; ++i) {
    header.totalControlPointCount += curves_[i].get_q().size();
  }
  size_t writesize;
  // Write header
  writesize = fwrite(&header, sizeof(header), 1, fp);
  if (writesize != 1) {
    std::cerr << "Cannot write BCC file header: " << filename << ".\n"
              << std::endl;
    throw std::runtime_error("Unable to load BCC file.");
  }

  Eigen::Matrix3Xf point_buffer;
  // Write curves
  for (int i = 0; i < header.curveCount; ++i) {
    // Write curve control point count
    int32_t curveControlPointCount = curves_[i].get_q().size();
    if (curves_[i].is_looped())
      curveControlPointCount = -curveControlPointCount;
    writesize = fwrite(&curveControlPointCount, sizeof(int32_t), 1, fp);
    if (writesize != 1) {
      std::cerr << "Cannot write curve control point count: " << filename
                << ", curve number " << i << ".\n"
                << std::endl;
      throw std::runtime_error("Unable to load BCC file.");
    }
    if (curves_[i].get_q().size() != 0) {
      Eigen::Map<const Eigen::Matrix3Xd> uncasted_points(
          curves_[i].get_q()[0].data(), 3, curves_[i].get_q().size());
      point_buffer = uncasted_points.cast<float>();
      writesize = fwrite(point_buffer.data(), sizeof(float),
                         3 * curves_[i].get_q().size(), fp);
      if (writesize != 3 * curves_[i].get_q().size()) {
        std::cerr << "Cannot write all points: " << filename
                  << ", curve number " << i << ", wrote " << writesize << " of "
                  << curves_[i].get_q().size() << ".\n"
                  << std::endl;
        throw std::runtime_error("Unable to load BCC file.");
      }
    }
  }

  fclose(fp);
}

} // namespace VerifyCurves