SYSTEM = M_species.o M_atom_functions.o M_kpoints.o M_atomPP_functions_KB.o M_configuraciones_McWEDA.o make_cells.o\
M_integrals_3c.o M_neighbors_PP.o M_atomPP_ion_functions_KB.o M_integrals_2c.o M_neighbors.o M_neighbors_vdW.o M_rotations.o

all: libthd_SYSTEM.a

.PHONY: clean

clean:
	rm -rf *.a *.o *.mod *genmod.f90 *~

libthd_SYSTEM.a: $(SYSTEM)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

#libthd_SYSTEM_CONFIG.a: $(SYSTEM_CONFIG)
#	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

%.o: %.f90
	${F90} ${FFLAGS} -c $< -I../include
