#!/usr/bin/env sh
set -ex

# drawbicyc in RelWithDebInfo mode
cd drawbicyc
mkdir build
cd build
cmake ..
VERBOSE=1 make -j4
cd ../..

set +ex # see https://github.com/travis-ci/travis-ci/issues/6522
