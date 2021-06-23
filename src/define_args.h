/**
 * Copyright (C) 2021 Ante Qu <antequ@cs.stanford.edu>.
 */

#pragma once

#include <gflags/gflags.h>

DEFINE_string(
    input, "data/glove.bcc",
    "The filepath to the input model, which must be in the BCC (Binary Curve "
    "Collection) file format.");

DEFINE_string(
    output, "results/out.txt",
    "The filepath for the desired output certificate. The first line of the "
    "output will be the number of curves. The remaining rows are a list of "
    "triplets, one comma-deliminated triplet per row. Each triplet contains a "
    "pair of curve indices, followed by the linking number between them. Pairs "
    "with zero linkage are omitted.");

DEFINE_string(
    comparewith, "[none]",
    "To verify against an existing topology certificate, set this option to "
    "the filepath to the input certificate. Otherwise, keep this set to "
    "'[none]'.");

DEFINE_string(
    diffpairs, "[none]",
    "If 'comparewith' is set, you can set this to an output filepath for a "
    "list of inconsistent loop pairs. The output rows are a list of triplets, "
    "one comma-deliminated triplet per row. Each triplet consists of a pair of "
    "curve indices, followed by a value that represents the difference between "
    "the two certificates for this pair. Keep this set to '[none]' to omit "
    "this output.");

DEFINE_string(
    diffcurves, "[none]",
    "If 'comparewith' is set, you can set this to an output filepath for a "
    "list of loops that participate in inconsistent loop pairs one loop index "
    "per row. Keep this set to '[none]' to omit this output.");

DEFINE_string(method, "barneshut",
              "Specifies the method for computing linking numbers. Possible "
              "values are 'barneshut' and 'directsum', and this is case "
              "sensitive.");

DEFINE_double(barnes_hut_init_beta, 2.0,
              "Advanced: sets the initial beta parameter for Barnes–Hut.");

DEFINE_double(barnes_hut_beta_limit, 10.0,
              "Advanced: sets the maximum beta parameter for Barnes–Hut.");
