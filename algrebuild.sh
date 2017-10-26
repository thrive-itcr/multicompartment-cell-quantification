#!/bin/bash
# Copyright (c) General Electric Company, 2017.  All rights reserved.

echo ""
echo "Removing algorithm image."
docker rmi rt106/multicompartment-cell-quant

echo ""
echo "Building algorithm image."
docker build -t rt106/multicompartment-cell-quant --build-arg http_proxy=$http_proxy --build-arg https_proxy=$https_proxy .

echo ""
docker images
