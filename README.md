# Fast Linking Numbers

<p align="center"><img src="https://graphics.stanford.edu/papers/fastlinkingnumbers/assets/RepImage.jpg" alt="Fast Linking Numbers" width="480"/></p>

This tool, called `verifycurves`, takes input models that consist of closed-loop curves, and outputs a topology certificate as a `.txt` file. It is based on [Qu and James 2021]. Among the different linking-number computation methods described in [Qu and James 2021], this tool supports using Barnes–Hut or direct summation to evaluate the Gauss linking integral.

Project Webpage: https://graphics.stanford.edu/papers/fastlinkingnumbers/

Supplemental: Derivation of the Two-Tree Barnes–Hut expansion terms used in this implementation [PDF](https://graphics.stanford.edu/papers/fastlinkingnumbers/assets/Qu2021BarnesHutExpandedTermsReference.pdf)

Original Code: https://github.com/antequ/fastlinkingnumbers

This branch is solely built for the purposes of reproducing Tables 1 and 2 in the paper, for the [Graphics Replicability Stamp Initiative](http://www.replicabilitystamp.org/).

## Instructions

We ran this on a Ubuntu 18.04 on a single 18-core Intel Xeon E5-3697 v4 @ 2.30 GHz with 64GB of RAM.

To perform the initial unzipping, you'll need to run the following commands (also in `unzip_and_run_scripts.sh`):

```
sudo apt-get install unzip -y
unzip fastlinkingnumbers.zip
cd fastlinkingnumbers
```

Afterwards, you can just run the following command to finish the set up and generate the tables.
```
. replicability.sh
```

The results will be in `table1.csv` and `table2.csv`. We've provided a copy from our machine, in `antequ_table1.csv` and `antequ_table2.csv`.

The script for generating these tables is in `src/replicability_script.cpp`.

## A Few Notes

1. We omitted Counting Crossings, FMM, and the GPU (DSG and BHG) columns in Table 1 because their licensing is more complex to release. We just provide Direct Summation (DS) and Barnes-Hut (BH).

2. There are small discrepencies in P and N_D for many inputs because some of them were originally measured using a single-precision implementation, while we released a double-precision implementation. Furthermore, some yarn examples differ a little more because they had been converted into the public-release BCC file-format in a process that converted between double and single precision. Despite these differences, the final linking numbers match what we had in the original paper, and the timing and error values should be on the same order as those in the paper. This is because the linking number is robust to variations in PLS and discretization, as long as both are performed correctly.

## Citation

Ante Qu and Doug L. James. 2021. Fast Linking Numbers for Topology Verification of Loopy Structures. *ACM Trans. Graph.* 40, 4, Article 106 (August 2021), 19 pages. https://doi.org/10.1145/3450626.3459778 

## License

The [MIT license](https://mit-license.org/) applies to all files in this codebase. In summary, you can reuse, modify, or distribute this, provided that you give attribution to us, at the minimum by copying the `LICENSE` file.

In addition, the [MPL2 license](https://www.mozilla.org/en-US/MPL/2.0/) applies to the following three files derived from the [Eigen](https://gitlab.com/libeigen/eigen) codebase:
* `src/BVAlgorithms.h`
* `src/BVH.h`
* `src/ModifiedKdBVH.h`

If you modify them, you must maintain the license text in the headers of these files. If you have questions, reference their FAQ [here](https://www.mozilla.org/en-US/MPL/2.0/FAQ/).

## Additional Attribution

The knitted yarn models are modified from Cem Yuksel's Yarn-level Cloth Models located at http://www.cemyuksel.com/research/yarnmodels/.
