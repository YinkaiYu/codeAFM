FC = mpiifort
FLAGS = -c -O3
LF = -qmkl
HOME = ../../../Lib_90_new
LIBS = $(HOME)/Modules/modules_90.a \
       $(HOME)/MyEis/libeis.a \
       $(HOME)/MyNag/libnag.a \
       $(HOME)/MyLin/liblin.a \
       $(HOME)/Ran/libran.a

all:
	cp $(HOME)/Modules/*.mod . ;\
	(make -f Compile  FC="$(FC)" `$(FC) -showme:compile`  LF="$(LF)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)")

clean: 
	(make -f Compile  clean );\
	rm *.mod