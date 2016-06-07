FDATA= M_Fdata_2c.o M_Fdata_3c.o
FDATA_DOGS = M_Fdata_1c_McWEDA_DOGS.o
FDATA_HARRIS = M_Fdata_1c_McWEDA_Harris.o

all: libthd_FDATA.a libthd_FDATA_Dogs.a libthd_FDATA_Harris.a

.PHONY: clean

clean:
	rm -rf *.a *.o *.mod *genmod.f90 *~

libthd_FDATA.a: $(FDATA)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

libthd_FDATA_Dogs.a: $(FDATA_DOGS)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

libthd_FDATA_Harris.a: $(FDATA_Harris)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

%.o: %.f90
	${F90} ${FFLAGS} -c $< -I../include
