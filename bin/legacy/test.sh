#! /bin/bash
#
gfortran -c stochastic_rk.f
gfortran -c -cpp stochastic_rk_test.f
gfortran -o stochastic_rk_test stochastic_rk_test.o stochastic_rk.o
rm stochastic_rk_test.o stochastic_rk.o
./stochastic_rk_test
rm stochastic_rk_test
