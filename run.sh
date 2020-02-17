#!/bin/bash

for i in {1..24}; do
    mkdir $i
    mkdir $i/DataFiles
    rsync ../Katana_Files*.m $i/
    rsync ../Katana_Files/DataFiles/*.nc $i/DataFiles
    rsync ../Katana_Files/Matlab_Startup $i/Matlab_Startup
done

