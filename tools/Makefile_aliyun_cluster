# the dir.where netcdf.mod exist
NETCDF  = /usr/lib64/gfortran/modules
MPIRUN  = /opt/openmpi/1.10.7/bin/mpirun

FC   =  /opt/openmpi/1.10.7/bin/mpifort

FF = -O2 -g -Wall -L/usr/lib64 -lnetcdf -lnetcdff -I$(NETCDF) 

EXE  = main
OBJS = mod_kind.o \
	mod_type.o \
	mod_param.o mod_con.o \
	mod_arrays.o \
	mod_io.o \
	mod_mympi.o \
	mod_debug.o \
	mod_den.o \
	mod_op.o \
	mod_int.o \
	main.o

.PHONY: all compile run clean

# need this 'empty' line to avoid some problem when 'make'
.SUFFIXES:

.SUFFIXES: .f90 .o

all: clean compile run view

compile: $(OBJS)
	$(FC) $(FF) $(OBJS) -o $(EXE)

.f90.o:
	$(FC) $(FF) -c $<

run:
	$(MPIRUN) -n 4 ./$(EXE)

view:
	ncview output/test.nc

clean:
	rm -f *.mod *.o *.out $(EXE)
