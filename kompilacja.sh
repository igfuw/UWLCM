cd build_LT
rm -rf *
cmake .. -Dlibcloudph++_DIR=/home/pzmij/biblioteki/local_folder/coal/share/libcloudph++ -Dlibmpdata++_DIR=/home/pzmij/biblioteki/local_folder/RWDI/share/libmpdata++ -DCMAKE_INSTALL_PREFIX=/home/pzmij/biblioteki/local_folder/RWDI -DCMAKE_BUILD_TYPE=RelWithDebInfo
make -j 20
make install -j20
