#!env bash
set -o pipefail

R -e "pak::pkg_install('bhklab/CoreGx')"
R -e "pak::pkg_install('bhklab/PharmacoGx')"