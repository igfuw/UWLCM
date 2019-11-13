Use Singularity: https://sylabs.io/singularity/ to create image for UWLCM

Build the image remotely here: https://cloud.sylabs.io/builder

```bash
$ singularity build --remote sng_ubuntu_18_04_cuda_10_0.sif sng_ubuntu_18_04_cuda_10_0
```

The image works with CUDA Driver Version: 410.79 and CUDA Version: 10.0

To run the image with GPU support use --nv flag:
```bash
$ singularity shell --nv sng_ubuntu_18_04_cuda_10_0.sif
```

It may be necessary to run the image in a clean environment: 
```bash
$ env -i singularity shell --nv sng_ubuntu_18_04_cuda_10_0.sif
```

Follow the libmpdata and libcloud readme on how to compile them and run tests.

Install both libraries in a local_folder:
```bash
$ cmake .. -DCMAKE_INSTALL_PREFIX=/local_folder/
$ make install
```

When compiling UWLCM tell cmake where to find the two locally installed libraries:
```bash
$ cmake .. -Dlibmpdata++_DIR=/local_folder/usr/local/share/libmpdata++ -Dlibcloudph++_DIR=/local_folder/usr/local/lib/cmake/libcloudph++
$ make
```

Prepare to wait as UWLCM and libcloudph++ take long time to compile.
After that youre done.
