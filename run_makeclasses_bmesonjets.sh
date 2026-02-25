#!/bin/bash

root -q macro_makeclasses_bmesonjets.cpp
mv T*.C T*.h ./include

echo "NOTE : Remember to set the mother and input folders in the directories.h file in the include folder!!"
