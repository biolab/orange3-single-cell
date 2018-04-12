#!/usr/bin/env bash

INSTALLER="$( cd "$(dirname "$0")/.." ; pwd -P )"
BUILD="$INSTALLER/windows/build"
DIST="$( cd "$(dirname "$0")/../../dist" ; pwd -P )"

# Build updated conda package
mkdir -p "$BUILD/conda"
rm -rf "$BUILD/conda"/*
conda build "$INSTALLER/conda" \
  --output-folder "$BUILD/conda"

# Build new conda package
NEW_SPEC="$BUILD/conda/conda-spec.txt"
sed '/^file/ d' \
  "$INSTALLER/windows/specs/conda-spec.txt" \
  > $NEW_SPEC 
CONDA_PACKAGE=$( find "$BUILD/conda" \
  -name "orange3-single-cell*" \
  -exec echo "file://{}" \; )
echo "$CONDA_PACKAGE" >> $NEW_SPEC

# Build new installer
./build-conda-installer.sh \
  --env-spec $NEW_SPEC \
  --online=no \
  --dist-dir $DIST

# Sign the installer
VERSION=$( echo $CONDA_PACKAGE | sed -n 's/.*-\([0-9.]*\)-.*/\1/p' )
signcode \
  -spc ~/Desktop/ulfri.spc \
  -v ~/Desktop/ulfri.pvk \
  -a sha1 \
  -t http://timestamp.verisign.com/scripts/timstamp.dll \
  -n scOrange \
  -i http://singlecell.biolab.si \
  "$DIST/scOrange-$VERSION-Miniconda-x86_64.exe"
