#!/bin/bash
R -e "library('devtools')" -e "devtools::build(binary = TRUE, args = c('--preclean'))"
R -e "install.packages(\"../rgoslin_2.0.1.1000_R_x86_64-pc-linux-gnu.tar.gz\", repos = NULL, type=\"source\")"
