#!/bin/bash

root -b -q "MakeCorrections.C(\"nominal\")"
root -b -q "MakeCorrections.C(\"jetid\")"
root -b -q "MakeCorrections.C(\"jesjer\")"
root -b -q "MakeCorrections.C(\"prior\")"
root -b -q "MakeCorrections.C(\"recseleff\")"
root -b -q "MakeCorrections.C(\"fitsignal\")"
