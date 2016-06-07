NAMD=

all: libthd_NAMD.a 

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
