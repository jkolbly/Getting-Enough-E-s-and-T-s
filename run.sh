gfortran main.f90 -o tmp > out 2>&1
./tmp "$@"
rm tmp