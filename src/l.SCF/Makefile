SCF = M_charges_Lowdin.o Qmixer_Anderson_DOGS.o Qmixer_Anderson_Harris.o
# Qmixer_Anderson_grid.o

all: libthd_SCF.a

.PHONY: clean

clean:
	rm -rf *.a *.o *.mod *genmod.f90 *~

libthd_SCF.a: $(SCF)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

%.o: %.f90
	${F90} ${FFLAGS} -c $< -I../include
