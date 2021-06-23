#!/bin/bash


echo "Glove: (Mean error 0.00214443)"
name=glove
./verifycurves --input=data/${name}.bcc --output=results/${name}.txt > ${name}_stdout.txt

if cmp -s results/${name}.txt results/reference_certificates/${name}.txt; then
  echo "Certificate matches the reference for ${name}."
else
  echo "Certificate does not match the reference for ${name}."
fi

echo "Knit Tube Init: (Mean error 0.00366683)"
name=knittubeinit
./verifycurves --input=data/${name}.bcc --output=results/${name}.txt > ${name}_stdout.txt

if cmp -s results/${name}.txt results/reference_certificates/${name}.txt; then
  echo "Certificate matches the reference for ${name}."
else
  echo "Certificate does not match the reference for ${name}."
fi

echo "Knit Tube Broken: (Mean error 0.006343)"
name=knittubebroken
./verifycurves --input=data/${name}.bcc --output=results/${name}.txt > ${name}_stdout.txt

if cmp -s results/${name}.txt results/reference_certificates/${name}.txt; then
  echo "Certificate matches the reference for ${name}."
else
  echo "Certificate does not match the reference for ${name}."
fi

echo "Simplified Woundball (blenderobj): (Mean error 0.000668562)"
name=blenderobj
./obj2bcc --input data/${name}.obj --output data/${name}.bcc > ${name}_obj2bcc_stdout.txt
./verifycurves --input data/${name}.bcc --output results/${name}.txt > ${name}_stdout.txt

if cmp -s results/${name}.txt results/reference_certificates/${name}.txt; then
  echo "Certificate matches the reference for ${name}."
else
  echo "Certificate does not match the reference for ${name}."
fi

echo "Simplified Woundball (houdiniobj): (Mean error 0.00138268)"
name=houdiniobj
./obj2bcc --input data/${name}.obj --output data/${name}.bcc > ${name}_obj2bcc_stdout.txt
./verifycurves --input data/${name}.bcc --output results/${name}.txt > ${name}_stdout.txt

if cmp -s results/${name}.txt results/reference_certificates/${name}.txt; then
  echo "Certificate matches the reference for ${name}."
else
  echo "Certificate does not match the reference for ${name}."
fi
