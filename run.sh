#!/bin/bash
cd ..
make
cd result
rm *result*
./../dynearthsol3d ./core-complex-3d-meshopt.cfg
python ../2vtk.py -t result
