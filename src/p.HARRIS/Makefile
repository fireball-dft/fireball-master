HARRIS=M_bcna_3c_Harris.o M_Coulomb.o M_kinetic.o M_overlap.o M_vna_Harris.o M_rho_2c_McWEDA_Harris.o M_rho_3c_McWEDA_Harris.o M_rhoS_2c_McWEDA_Harris.o M_rhoS_3c_McWEDA_Harris.o M_vnl_KB.o
VXC_MCWEDA = M_vxc_McWEDA_Harris.o
VXC_HORSFIELD = M_vxc_Horsfield_Harris.o

all: libthd_HARRIS.a

.PHONY: clean

clean:
	rm -rf *.a *.o *.mod *genmod.o *~

libthd_HARRIS.a: $(HARRIS)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

libthd_VXC_MCWEDA.a: $(VXC_MCWEDA)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

libthd_VXC_HORSFIELD.a: $(VXC_HORSFIELD)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

%.o: %.f90
	${F90} ${FFLAGS} -c $< -I../include
