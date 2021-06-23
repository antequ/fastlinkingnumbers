# Fast Linking Numbers

<p align="center"><img src="https://graphics.stanford.edu/papers/fastlinkingnumbers/assets/RepImage.jpg" alt="Fast Linking Numbers" width="480"/></p>

This tool, called `verifycurves`, takes input models that consist of closed-loop curves, and outputs a topology certificate as a `.txt` file. It is based on [Qu and James 2021]. Among the different linking-number computation methods described in [Qu and James 2021], this tool supports using Barnes–Hut or direct summation to evaluate the Gauss linking integral.

Project Webpage: https://graphics.stanford.edu/papers/fastlinkingnumbers/

Supplemental: Derivation of the Two-Tree Barnes–Hut expansion terms used in this implementation [PDF](https://graphics.stanford.edu/papers/fastlinkingnumbers/assets/Qu2021BarnesHutExpandedTermsReference.pdf)

Code: https://github.com/antequ/fastlinkingnumbers

## Installation

For Ubuntu, we tested with Ubuntu 18.04, gcc v7.5.0, CMake v3.10.2, gflags v2.2.1, libomp v10.0.0, and Eigen v3.3.4. For Windows, we tested with Windows 10, Visual Studio Community 2019 v16.10.2, CMake v3.20.5, and gflags v2.2.2, for 64-bit release.

### Required packages:
* [CMake](https://cmake.org/)
* [Eigen 3](https://eigen.tuxfamily.org/)
* [gflags](https://gflags.github.io/gflags/)
* OpenMP 3.1

You can install them on Ubuntu 18.04 via
```
sudo apt-get install libgflags-dev libeigen3-dev libomp-10-dev build-essential cmake -y
```

### Building the executable verifycurves tool on Ubuntu

From the root source directory, you can build this like any other CMake project with these commands:
```
mkdir build
cd build
cmake ..
make -j
cd ..
ln -s build/verifycurves verifycurves
ln -s build/obj2bcc obj2bcc
```

To build only the `verifycurves` executable, use `make -j verifycurves`. Similarly you can use `make -j obj2bcc` or `make -j libverifycurves` to build the file format converter or the dynamic-link library.

### Running on Windows

Precompiled Binaries:

We've uploaded the compiled 64-bit windows binaries [here](https://drive.google.com/file/d/1nLfxRrf1zoVjY3bT1U-xhz6KByn7yd-S/view?usp=sharing). You may need to install the [Microsoft Visual C++ Redistributable for Visual Studio 2019](https://visualstudio.microsoft.com/downloads/#microsoft-visual-c-redistributable-for-visual-studio-2019) runtime libraries if you do not already have them.

Building from the source code:

To build from the source code, download and install [CMake](https://cmake.org/), [Eigen 3](https://eigen.tuxfamily.org/), [gflags](https://gflags.github.io/gflags/), and Microsoft Visual Studio 2019. Make a `build` subdirectory. In `cmake-gui`, set the project source directory and the build subdirectory. Click "Configure", making sure that it finds Eigen3 and gflags. Then click "Generate", then "Open Project". Build the solution, for Release x64, in Visual Studio.

### Shared dynamic-link library

Alternately, you can also dynamically link this project against your code. To do so on Ubuntu, run the following lines after the above:
```
cd build
sudo make install
cd ..
```

Then simply include `includes/curve.h`, `includes/model.h` and `includes/linking_number_certificate.h`, as desired, into your project. Note that all of our projects are standalone and do not depend on this shared library.


## Usage

### Computing a topology certificate for an input model
Run `verifycurves` on any input model in the [BCC file format](#curve-input-file-format) to generate a certificate. The input models must contain only closed loops. The certificate will be a sparse matrix in triplet form, in plaintext (`.txt`). The first row is the number of curves in the model, and in every subsequent row, the first two coordinates are the curve indices, and the third is the linking number between them. Curve indices are zero-based. Example input files are provided in the `data` folder. Here is an example command:

```
./verifycurves --input data/glove.bcc --output results/glove.txt
```

This should be the output:

`glove.txt`:
```
70
```

This output indicates that there are 70 curves and no pair of curves are linked.

You can also set `--method barneshut` or `--method directsum` to use Barnes–Hut or Direct Summation to compute the linking numbers. This tool by default uses Barnes–Hut. Here is another example command that tells the tool to use Direct Summation:

```
./verifycurves --input data/knittubebroken.bcc --output results/knittubebroken.txt --method directsum
```

This should be the output:

`knittubebroken.txt`:
```
39
5,6,-2
```

This output indicates that there are 39 curves, and curves 5 and 6 have a linking number of -2 between them, while all other curve pairs have a linking number of 0.


### Comparing two topology certificates
As these sparse matrix representations are unique, one way to check the certificates is just using the bash `diff` command on the two output certificate files.
```
./verifycurves --input data/knittubeinit.bcc --output results/knittubeinit.txt
./verifycurves --input data/knittubebroken.bcc --output results/knittubebroken.txt
diff results/knittubeinit.txt results/knittubebroken.txt
```

### Verifying an input model against an existing topology certificate

Alternately, you can also just run `verifycurves` on an input model with the argument `--comparewith [CERTIFICATE_FILENAME]` to compare the input model against the existing certifcate in `[CERTIFICATE_FILENAME]`. You can also have it output a list of inconsistent loop pairs via `--diffpairs` and a list of loops that participate in these pairs via `--diffcurves`.

```
./verifycurves --input data/knittubebroken.bcc --comparewith results/knittubeinit.txt --diffpairs results/knittube.dp.txt --diffcurves results/knittube.dc.txt
```

These should be the output:

`knittube.dp.txt`:
```
5,6,-2
```

`knittube.dc.txt`:
```
5
6
```
These outputs indicate, in this example, that the link between curves 5 and 6 is incorrect, with a difference of -2 between the before and after links.

### Curve input file format

`verifycurves` takes the simple BCC (Binary Curve Collection) file format as input. See the file format specification from Cem Yuksel's Yarn-level Models webpage [here](http://www.cemyuksel.com/research/yarnmodels/) and an example loading tutorial [here](http://www.cemyuksel.com/cyCodeBase/soln/using_bcc_files.html). In addition to uniform Catmull–Rom splines, we added three more curve formats, so Bytes 4–5 in the header can be any of the following:
* `C0` - [Catmull–Rom curves with uniform parameterization](http://www.cemyuksel.com/research/catmullrom_param)
* `BS` - Cubic [B-Splines](https://en.wikipedia.org/wiki/B-spline)
* `BZ` - Cubic [Bézier curves](https://en.wikipedia.org/wiki/B%C3%A9zier_curve)
* `PL` - [Polylines](https://en.wikipedia.org/wiki/Polygonal_chain)

The `data` folder has example BCC data files. You can find more input loopy structures from Cem Yuksel's [Yarn-Level Models](http://www.cemyuksel.com/research/yarnmodels/) (the first 8 are looped knit structures that we can handle) and our extended input dataset (described [below](#extended-input-dataset)).

### Converting from Blender and Houdini Obj files

Blender and Houdini curves can be passed to our tool by exporting them as lines in the OBJ file format, and then using `obj2bcc` to convert it into BCC as closed curves.

Run `obj2bcc` with the `--input [INPUT]` argument for input, `--segment_type [TYPE]` to specify the segment type (we do not use the `cstype` field in the OBJ format), and `--output [OUTPUT]` for output:

```
./obj2bcc --input data/houdiniobj.obj --segment_type PL --output houdiniobj.bcc
```

Within the OBJ files, control points or polyline vertices must be represented as "l" (lines) rather than "curv" (curves) or "f" (faces). While most Houdini OBJ files work out of the box, when exporting from Blender, we require these line segments to be in order. One way to ensure this from Blender is by [converting the object](https://docs.blender.org/manual/en/latest/scene_layout/object/editing/convert.html) first to "Grease Pencil", and then converting it again to "Polygonal Curve." Then when you export as OBJ, the line segments will be in order. We've provided an example Blender obj and an example Houdini obj in the `data` folder, as `blenderobj.obj` and `houdiniobj.obj`, respectively; these both correspond to the simpler Woundball used for the Figure 6 illustration.

### Testing

We've included a few example data files in `data/` along with their reference certificates in `results/reference_certificates`. You can run `. test.sh` to generate your own output certificates for these files to compare against the reference certificates. Our extended dataset also has reference certificates included.

Sadly, we do not release unit or other tests in this code release. Modify this code at your own risk, and we suggest you create your own tests.

### Accuracy versus performance tradeoff

The initial settings we provided produce correct linking numbers for our entire dataset and most yarn simulation data we have encountered. However, for certain specific other curve collections, you may need to tune the parameters to increase the accuracy if you do not trust the results. Here are a few ways:

To increase accuracy at the cost of performance, you can manually increase the `--barnes_hut_init_beta` and/or `--barnes_hut_beta_limit` parameters from their current default values of 2.0 and 10.0. Increasing beta drastically reduces error of the first Barnes–Hut run, but increases its runtime. For example, try:
```
./verifycurves --input data/glove.bcc --output results/glove.txt --barnes_hut_init_beta 5.0 --barnes_hut_beta_limit 10.0
```

The most accurate, yet slowest, method is to just use direct summation:
```
./verifycurves --input data/glove.bcc --output results/glove.txt --method directsum
```


## Extended input dataset

Example input data can be downloaded from the Fast Linking Numbers [project webpage](https://graphics.stanford.edu/papers/fastlinkingnumbers/). In particular, the download link is [here (ZIP, 980 MB)](https://drive.google.com/file/d/1tTSrzwP92xYxmVXVc6bx3GGEJ0srukp9/view?usp=sharing) for the dataset and [here (TXT)](https://graphics.stanford.edu/papers/fastlinkingnumbers/assets/dataset_readme.txt) for the dataset readme.

## Citation
Ante Qu and Doug L. James. 2021. Fast Linking Numbers for Topology Verification of Loopy Structures. *ACM Trans. Graph.* 40, 4, Article 106 (August 2021), 19 pages. https://doi.org/10.1145/3450626.3459778 

## License
The [MIT license](https://mit-license.org/) applies to all files in this codebase. In summary, you can reuse, modify, or distribute this, provided that you give attribution to us, at the minimum by copying the `LICENSE` file.

In addition, the [MPL2 license](https://www.mozilla.org/en-US/MPL/2.0/) applies to the following three files derived from the [Eigen](https://gitlab.com/libeigen/eigen) codebase:
* `src/BVAlgorithms.h`
* `src/BVH.h`
* `src/ModifiedKdBVH.h`

If you modify them, you must maintain the license text in the headers of these files. If you have questions, reference their FAQ [here](https://www.mozilla.org/en-US/MPL/2.0/FAQ/).
