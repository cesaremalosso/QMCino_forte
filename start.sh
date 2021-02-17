gfortran -c moduli.f90 dmc_mw.f90 -O2 -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic  -fbacktrace
gfortran moduli.o dmc_mw.o -o go.exe
time ./go.exe
