#!/usr/bin/env sh
set -ex

# UWLCM in RelWithDebInfo mode without 3D and with 'abs' libmpdata++ option
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo -DUWLCM_DISABLE="2D_LGRNGN;2D_BLK_1M" -DMPDATA_OPTS="ABS"
VERBOSE=1 make
cd tests/moist_thermal
OMP_NUM_THREADS=6 make test || cat Testing/Temporary/LastTest.log / # "/" intentional! (just to make cat exit with an error code)
cd ../../..

set +ex # see https://github.com/travis-ci/travis-ci/issues/6522
