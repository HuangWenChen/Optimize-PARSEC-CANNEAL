# Optimize-PARSEC-CANNEAL

## HOW TO BUILD

```sh
cd ./src
```
```sh 
make
``` 
or can assign some params (default THREAD=1 COMPILER=g++)
```sh 
make THREAD=4 COMPILER=clang++-5.0
```

## HOT TO WORK
```
cd ./src
./canneal NTHREADS NSWAPS TEMP NETLIST [NSTEPS]
```

## HOW TO TESTS

On Linux,  use this makefile target to test file.
```sh
make test
```
If you want to test other testing file, you can put your testing file into 
inputs(directory) and run target
