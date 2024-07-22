#! /bin/sh

cd ../
mkdir build
cd build

rm -r windows_solution_x64_dll
mkdir windows_solution_x64_dll && cd windows_solution_x64_dll

cmake ../.. -Ax64 -DCMAKE_SYSTEM_NAME=windows -DBUILD_SHARED_LIBS=ON
