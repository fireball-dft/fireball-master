DOGS = M_bcna_3c_DOGS.o  M_dipole_z.o  M_vna_DOGS.o  M_vxc_McWEDA_DOGS.o

all: libthd_DOGS.a

.PHONY: clean

clean:
	rm -rf *.a *.o *.mod *genmod.o *~

libthd_DOGS.a: $(DOGS)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

%.o: %.f90
	${F90} ${FFLAGS} -c $< -I../include
