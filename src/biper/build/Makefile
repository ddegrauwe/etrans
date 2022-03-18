VPATH=../external:../module:../interface:../programs:../xrd/module
objects = test.o parkind1.o yomhook.o etibihie.o espline_mod.o esmoothe_mod.o fpbipere.o horiz_field.o
compiler = pgf95
lib_dir =../interface

test : $(objects) 
	$(compiler) -o test $(objects) -L. -lxrd-pgi -lmpidummy

test.o : parkind1.F90 yomhook.F90 espline_mod.F90 esmoothe_mod.F90 etibihie.F90 fpbipere.F90 horiz_field.F90 test.F90
	$(compiler) -c -I../interface -I../xrd/include $^

espline_mod.o : espline_mod.F90 
	$(compiler) -c  $^ 

etibihie.o : espline_mod.F90 esmoothe_mod.F90 etibihie.F90 
	$(compiler) -I../xrd/include -c  $^ 

fpbipere.o : espline_mod.F90 esmoothe_mod.F90 fpbipere.F90 
	$(compiler)  -c $^ 

horiz_field.o : horiz_field.F90 
	$(compiler)  -c  $^ 

esmoothe_mod.o : esmoothe_mod.F90 
	$(compiler)  -c  $^ 


clean: 
	rm test $(objects) *.mod *.txt


