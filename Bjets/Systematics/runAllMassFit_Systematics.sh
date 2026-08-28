#!/bin/bash

cd ../

#RecSelEff
root -l -b -q "MassFit.C(1, 1, 1, 0, 5.0, 250.0)";
#FitModel
root -l -b -q "MassFit.C(1, 1, 0, 1, 5.0, 250.0)";


cd Systematics/
