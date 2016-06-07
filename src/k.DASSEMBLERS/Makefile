DASSEMBLERS_HARRIS = M_DFdata_2c.o M_DFdata_3c.o M_Drotations.o M_Drotations_PP.o M_Dassemble_2c_Harris.o M_Dassemble_PP_2c.o M_Dassemble_rho_McWEDA_Harris.o M_Dassemble_usr_McWEDA_Harris.o M_build_forces.o M_Dassemble_3c_Harris.o M_Dassemble_vxc_McWEDA.o M_Dassemble_ewald_Harris.o
DASSEMBLERS_DOGS = M_DFdata_2c.o M_DFdata_3c.o M_Drotations.o M_Drotations_PP.o M_Dassemble_2c_DOGS.o M_Dassemble_PP_2c.o M_Dassemble_rho_McWEDA_DOGS.o M_Dassemble_usr_McWEDA_DOGS.o M_build_forces.o M_Dassemble_3c_DOGS.o M_Dassemble_vxc_McWEDA.o M_Dassemble_ewald_Harris.o

all: libthd_DASSEMBLERS_HARRIS.a libthd_DASSEMBLERS_DOGS.a

.PHONY: clean

clean:
	rm -rf *.a *.o *.mod *genmod.f90 *~

libthd_DASSEMBLERS_Harris.a: $(DASSEMBLERS_HARRIS)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

libthd_DASSEMBLERS_Dogs.a: $(DASSEMBLERS_DOGS)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

%.o: %.f90
	${F90} ${FFLAGS} -c $< -I../include
