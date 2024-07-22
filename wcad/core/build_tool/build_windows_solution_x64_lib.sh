#! /bin/sh

cd ../
mkdir build
cd build

rm -r windows_solution_x64_lib
mkdir windows_solution_x64_lib && cd windows_solution_x64_lib

cmake ../.. -Ax64 -DCMAKE_SYSTEM_NAME=windows -DBUILD_SHARED_LIBS=OFF
