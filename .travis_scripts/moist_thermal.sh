#!/usr/bin/env sh
set -ex

# build UWLCM in RelWithDebInfo mode without 3D and with 'abs' libmpdata++ option
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo -DUWLCM_DISABLE="3D_LGRNGN;3D_BLK_1M;PIGGYBACKER;SGS" -DMPDATA_OPTS="ABS"
VERBOSE=1 make -j2
sudo make install
cd ../..

# install UWLCM_plotting
git clone --depth=1 git://github.com/igfuw/UWLCM_plotting.git
cd UWLCM_plotting/drawbicyc
mkdir build
cd build
cmake ..
sudo make install
cd ../../../

# build UWLCM tests
cd tests
mkdir build && cd build
cmake ..
VERBOSE=1 make -j2

# run the moist_thermal test
cd moist_thermal
OMP_NUM_THREADS=6 make test || cat Testing/Temporary/LastTest.log / # "/" intentional! (just to make cat exit with an error code)
cd ../../..

set +ex # see https://github.com/travis-ci/travis-ci/issues/6522
