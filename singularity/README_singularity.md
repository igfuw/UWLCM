Use Apptainer or Singularity

Built images are available at: https://zenodo.org/records/15591478 https://cloud.sylabs.io/library/pdziekan/

Alternatively, build you own image using one of the .def files in this folder:

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

Within the image shell, follow the libmpdata and libcloud readme on how to compile them and run tests, e.g.:
```bash
$ cmake .. -DCMAKE_INSTALL_PREFIX={INSTALL_DIR}
$ make -j4 install
```

When compiling UWLCM tell cmake where to find the two locally installed libraries:
```bash
$ cmake .. -Dlibmpdata++_DIR={INSTALL_DIR}/usr/local/share/libmpdata++ -Dlibcloudph++_DIR={INSTALL_DIR}/usr/local/lib/cmake/libcloudph++
$ make -j4 install
```

Running the model (using the binary uwlcm) also has to be done in the shell.

