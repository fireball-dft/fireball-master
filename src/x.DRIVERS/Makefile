BINS = begin.x begin-looprc.x create.x fireball.x lightning.x

all: ${BINS}

.PHONY: clean

clean:
	rm -rf *.a *.o *.mod *genmod.f90 *~

%.o: %.f90
	${F90} ${FFLAGS} -c $< -I../include
