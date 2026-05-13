Use Apptainer or Singularity

Built images are available at: https://zenodo.org/records/15591478 https://cloud.sylabs.io/library/pdziekan/

Alternatively, build you own image using one od the .def files in this folder

```bash
$ singularity build X.sif X.def
```

To run a shell within the image (the --nv flag is needed for GPUs to work):
```bash
$ singularity shell --nv X.sif
```

It may be necessary to run the image in a clean environment: 
```bash
$ env -i singularity shell --nv X.sif
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

