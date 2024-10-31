#!/bin/bash

find /home/connor/dev/geostan/tests -maxdepth 1 -type f -exec Rscript {} \; > test_results.out


