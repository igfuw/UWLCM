#!/usr/bin/env sh
set -ex

# UWLCM in RelWithDebInfo mode
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo -DUWLCM_DISABLE="SGS"
VERBOSE=1 make
cd tests/unit
OMP_NUM_THREADS=1 ctest -R api_test_iles || cat Testing/Temporary/LastTest.log / # "/" intentional! (just to make cat exit with an error code)
cd ../../..

set +ex # see https://github.com/travis-ci/travis-ci/issues/6522
