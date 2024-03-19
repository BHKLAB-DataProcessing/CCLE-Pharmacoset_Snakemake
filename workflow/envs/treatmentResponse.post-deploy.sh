#!env bash
set -o pipefail

# This will force install the newest development version of the package
R -e "pak::pkg_install('bhklab/CoreGx')"
R -e "pak::pkg_install('bhklab/PharmacoGx')"