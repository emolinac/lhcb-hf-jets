#!/bin/bash

root -b -q "MassFit.C(0,\"nominal\")"
root -b -q "MassFit.C(0,\"recseleff\")"
root -b -q "MassFit.C(0,\"fitsignal\")"

root -b -q "MassFit.C(1,\"nominal\")"
root -b -q "MassFit.C(1,\"recseleff\")"
root -b -q "MassFit.C(1,\"fitsignal\")"