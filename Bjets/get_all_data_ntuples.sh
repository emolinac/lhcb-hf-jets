#!/bin/bash

root -b -q "MakeVarTreeData2016MU.C(0)"
root -b -q "MakeVarTreeData2016MD.C(0)"
root -b -q "MakeVarTreeData2017MU.C(0)"
root -b -q "MakeVarTreeData2017MD.C(0)"
root -b -q "MakeVarTreeData2018MU.C(0)"
root -b -q "MakeVarTreeData2018MD.C(0)"

root -b -q "MakeVarTreeData2016MU.C(1)"
root -b -q "MakeVarTreeData2016MD.C(1)"
root -b -q "MakeVarTreeData2017MU.C(1)"
root -b -q "MakeVarTreeData2017MD.C(1)"
root -b -q "MakeVarTreeData2018MU.C(1)"
root -b -q "MakeVarTreeData2018MD.C(1)"

cd ../output-files/

hadd ntuple_bjets_data.root ntuple_bjets_data_201*_M*_nominal.root

hadd ntuple_bjets_data_jetid.root ntuple_bjets_data_201*_M*_jetid.root

cd ../Bjets/