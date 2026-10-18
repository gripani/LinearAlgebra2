FC = gfortran
FFLAGS = -std=f2008 -fPIC -O3 -Wall -W -Werror -Wno-compare-reals
MODDIR = bin/mod

.PHONY: all clean

all: lib/linearalgebra_wrapper.so

$(MODDIR):
	mkdir -p $@

bin:
	mkdir -p $@

lib:
	mkdir -p $@

bin/linearalgebra.o: src/linearalgebra.f90 | bin $(MODDIR)
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

bin/linearalgebra_wrapper.o: src/linearalgebra_wrapper.f90 bin/linearalgebra.o | bin $(MODDIR)
	$(FC) $(FFLAGS) -I$(MODDIR) -J$(MODDIR) -c $< -o $@

lib/linearalgebra_wrapper.so: bin/linearalgebra.o bin/linearalgebra_wrapper.o | lib
	$(FC) $(FFLAGS) -shared $^ -o $@

clean:
	rm -f bin/*.o bin/mod/*.mod lib/linearalgebra_wrapper.so
