THEORY=M_Coulomb.f90  M_dipole_z.f90  M_kinetic.f90  M_overlap.f90 M_vxc_McWEDA_DOGS.f90 
VXC_MCWEDA = M_vxc_McWEDA_Harris.o
VXC_HORSFIELD = M_vxc_Horsfield_Harris.o

all: libthd_THEORY.a libthd_VXC_MCWEDA.a libthd_VXC_HORSFIELD.a

.PHONY: clean

clean:
	rm -rf *.a *.o *.mod *genmod.o *~

libthd_THEORY.a: $(THEORY)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

libthd_VXC_MCWEDA.a: $(VXC_MCWEDA)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

libthd_VXC_HORSFIELD.a: $(VXC_HORSFIELD)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

%.o: %.f90
	${F90} ${FFLAGS} -c $< -I../include
