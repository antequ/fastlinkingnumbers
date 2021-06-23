/**
 * Copyright (C) 2021 Ante Qu <antequ@cs.stanford.edu>.
 * Copyright (C) 2018 Jonathan Leaf.  All rights reserved.
 * Originally "sim_model.h" by Jonathan Leaf. Modified by Ante Qu.
 */

#pragma once
#include "curve.h"
#include <vector>

namespace VerifyCurves {
/**
 * A Model stores a collection of curves. This class has methods to import from
 * and export to the Binary Curve Collection (BCC) file format.
 *
 * This class also has a method, ImportFromObjFile, to assist in important and
 * converting from certain OBJ files. It can handle curve data from Blender or
 * Houdini, exported as OBJs, by reading vertices ("v") and lines ("l"). To
 * handle Blender OBJs, it will attach adjacent lines ("l") together into the
 * same curve if their last and first vertices respectively match. The caller
 * must specify the desired segment type for the curves. For best results, we
 * recommend working with the BCC file format.
 */
class Model {
public:
  /**
   * The default constructor generates an empty list of curves.
   */
  Model();

  /**
   * This method imports a curve collection from a BCC file and stores it in
   * this Model, in curves_.
   *
   * @param [in] filename the filepath of the BCC file to import from.
   *
   * @post curves_ is updated with the contents from the BCC file.
   */
  void ImportFromFile(const std::string &filename);

  /**
   * This method exports the curve collection in this Model to a BCC file.
   *
   * @param [in] filename the filepath of the BCC file to save to.
   */
  void ExportToFile(const std::string &filename) const;

  /**
   * These two methods get the list of curves, either as a modifiable or const
   * reference.
   */
  std::vector<Curve> &GetCurves() { return curves_; };
  const std::vector<Curve> &GetCurves() const { return curves_; };

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
  void ImportFromObjFile(const std::string &filename, SegmentType segment_type);

private:
  // These methods enable importing from and exporting to BCC files.
  void ImportBccFile(const std::string &filename);
  void ExportBccFile(const std::string &filename) const;

  // Curve data is stored in curves_.
  std::vector<Curve> curves_;
};

} // namespace VerifyCurves