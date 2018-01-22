NETCDF  = /opt/local
MPIRUN  = mpirun

FC   =  mpiifort

FF = -stack_temps -safe_cray_ptr -ftz -assume byterecl -fp-model precise -O2 -i4 -r8 -I$(NETCDF)/include -L$(NETCDF)/lib -lnetcdf -lnetcdff

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
	$(MPIRUN) -n 32 ./$(EXE)

view:
	ncview output/test.nc

clean:
	rm -f *.mod *.o *.out $(EXE)
