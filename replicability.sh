#!/bin/bash

echo "This code is from the replicability branch of the code release:"
echo "https://github.com/antequ/fastlinkingnumbers/tree/replicability,"
echo "and the dataset directory comes directly from the released dataset:"
echo ""

echo "Installing prerequisites (needs sudo access) ..."
sudo apt-get update
sudo apt-get install build-essential cmake libgflags-dev libeigen3-dev -y

echo "Building project ..."
mkdir build
cd build
cmake ..
make -j replicability
cd ..

echo "Running the replicability executable."
echo "This was built from the script in src/replicability_script.cpp."
echo "Tables 1 and 2 from the paper will be saved in table1.csv and table2.csv."
echo "The standard output will be piped to replicability_out.txt."
echo ""
echo "If this takes too long, feel free to hit CTRL+C. I've computed the results"
echo "on my machine and saved them as antequ_table1.csv and antequ_table2.csv"
echo "for you."
echo ""
echo "This takes about 15 minutes on a single 18-core Intel Xeon E5-3697 v4."
./build/replicability > replicability_out.txt

echo ""
echo "Table 1 of Results You Just Computed:"
cat table1.csv
echo ""
echo "Table 2 of Results You Just Computed:"
cat table2.csv
echo ""

echo ""
echo "Table 1 of Results Ante's Machine Computed:"
cat antequ_table1.csv
echo ""
echo "Table 2 of Results Ante's Machine Computed:"
cat antequ_table2.csv
echo ""


echo "A few notes:"
echo "1. We omitted Counting Crossings, FMM, and the GPU (DSG and BHG) columns"
echo "   in Table 1 because their licensing is more complex to release. We just"
echo "   provide Direct Summation (DS) and Barnes-Hut (BH)."
echo "2. There are small discrepencies in P and N_D for many inputs because some"
echo "   of them were originally measured using a single-precision"
echo "   implementation, while we released a double-precision implementation."
echo "   Furthermore, some yarn examples differ a little more because they had"
echo "   been converted into the public-release BCC file-format in a process"
echo "   that converted between double and single precision. Despite these"
echo "   differences, the final linking numbers match what we had in the"
echo "   original paper, and the timing and error values should be on the same"
echo "   order as those in the paper. This is because the linking number is"
echo "   robust to variations in PLS and discretization, as long as both are"
echo "   performed correctly."

