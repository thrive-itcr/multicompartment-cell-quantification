# Multi-Compartment Cell Quantification

__Copyright (c) General Electric Company, 2017.  All rights reserved.__


To build the Docker container:
```sh
$ docker build -t thrive/multi-compartment-cell-quantification --build-arg http_proxy=$http_proxy --build-arg https_proxy=$https_proxy --build-arg no_proxy=$no_proxy .
```
Run example:
python CellQuantificationMultiMarkers.py -Nucsegmask test/NucSeg.tif -Cellsegmask test/CellSeg.tif -inbioim test/test2_dapi.tif test/test2_pck26.tif -bioname dapi pck26 -outname test/out.csv
