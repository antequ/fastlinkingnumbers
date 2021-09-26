#!/bin/bash

# Run this from outside fastlinkingnumbers.zip
# on a clean Ubuntu build.

sudo apt-get install unzip -y
unzip fastlinkingnumbers.zip
cd fastlinkingnumbers

. replicability.sh
