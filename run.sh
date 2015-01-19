#!/bin/bash -f 

while [ $# -ne 1 ]; do 
      echo "Usage : $0 CASE "
      echo "   CASE  : evpprecondpcg  " 
      exit 0
done

CASE=$1

rm ./$CASE.exe
ifort -g -traceback -shared-intel -mcmodel=large -check  noarg_temp_created -convert big_endian -C ./$CASE.f90  -o ./$CASE.exe
chmod +x $CASE.exe
./$CASE.exe
