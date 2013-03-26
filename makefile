
OBJ1 = baseflow.f90 

code : $(OBJ1)
	mpif90 $(OBJ1) -o prg

clean :
	rm *bin *mod *~ coef fs *dat

.f.o:
	mpif90 -c -03 $<
