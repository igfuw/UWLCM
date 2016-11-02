#!/bin/sh
#BOOST_ROOT=/home/plgdziekan/local_installs/boost_1_60_0 
#PKG_CONFIG_PATH=/home/plgdziekan/local_installs/blitz-0.10 ~/local_installs/cmake-3.5.0-Linux-x86_64/bin/cmake .. -DCMAKE_CXX_COMPILER=$HOME/local_installs/gcc-4.9.3-install/bin/g++ -DCMAKE_BUILD_TYPE=Release
#PKG_CONFIG_PATH=/home/plgdziekan/local_installs/blitz-0.10 ~/local_installs/cmake-3.5.0-Linux-x86_64/bin/cmake .. -DCMAKE_PREFIX_PATH=/home/plgdziekan/local_installs/libcloudph++/usr/local -DCMAKE_CXX_COMPILER=$HOME/local_installs/gcc-4.9.3-install/bin/g++ -DCMAKE_BUILD_TYPE=Release
PKG_CONFIG_PATH=/home/plgdziekan/local_installs/blitz-0.10 ~/local_installs/cmake-3.5.0-Linux-x86_64/bin/cmake .. -DCMAKE_PREFIX_PATH="/home/plgdziekan/code/libmpdataxx;/home/plgdziekan/local_installs/libcloudph++/usr/local" -DCMAKE_CXX_COMPILER=$HOME/local_installs/gcc-4.9.3-install/bin/g++ -DCMAKE_BUILD_TYPE=Release

#-DCMAKE_LIBRARY_PATH=/home/plgdziekan/code/libcloudphxx/build/src
