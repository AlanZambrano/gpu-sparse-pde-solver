TAMANIO = 180

# ****descomentame para usar Magma***
MAGMAONOFF = -DMAGMA


## Paths ##
PATH_EIGEN = /state/partition1/apps/spack/opt/spack/linux-centos7-x86_64/gcc-11.1.0/eigen-3.4.0-xaed56heqf6sisylasau4wt7dvmbgp2r/include/eigen3/

## MAGMA ###
LIBSMAGMA=/share/apps/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/magma-2.7.1-lrzhx5oz5pjqikja2jz6czd3m47r5fzj/lib/

INCLUDEMAGMA=/share/apps/spack/opt/spack/linux-centos7-x86_64/gcc-8.2.0/magma-2.7.1-lrzhx5oz5pjqikja2jz6czd3m47r5fzj/include/

MAGMAFLAG = -I $(INCLUDEMAGMA) -L $(LIBSMAGMA) -lmagma -lmagma_sparse -lcudart -lcusparse -lcusolver -lcublas -std=c++11 -DHAVE_CUBLAS $(MAGMAONOFF)


########################
## Compiler and Flags ##

GPROF   = -pg

## para el cpu-prof descomenta el primero y comenta el segundo CC
#CC     = g++ $(GPROF)
CC     = nvcc -Wno-deprecated-gpu-targets


CFLAGS = $(MAGMAFLAG) -I $(PATH_EIGEN) -DTAMANIO=$(TAMANIO)
###########################################################


OBJS = main.o funciones.o estado.o

main: $(OBJS)
	$(CC) $(CFLAGS) $^ -o ../bin/main_$(TAMANIO)

main.o:	main.cpp ../include/funciones.h ../include/estado.h 
	$(CC) $(CFLAGS) -c $< -o $@

estado.o: estado.cpp estado_operadores.cpp estado_ensamblaje.cpp ../include/*.h
	$(CC) $(CFLAGS) -c $< -o $@

funciones.o: funciones.cpp ../include/parametros.h
	$(CC) $(CFLAGS) -c $< -o $@


# MAGMA con Eigen ####
magmatest: magmatest.o
	$(CC) $(MAGMAFLAG) -I $(PATH_EIGEN) $^ -o magmatest	

magmatest.o: magmatest.cu
	$(CC) $(MAGMAFLAG) -I $(PATH_EIGEN) -c $< -o $@
############

cusparse:
	$(CC) $(CFLAGS) cusparse_chatgpt.cu -o a.out

test:
	$(CC) $(CFLAGS) test.cu -o a.out

run:
	../bin/main

submit:
	qsub jobGPU; watch qstat

clean:
	rm *.o    #oveja* *.dat pred-p* ../bin/mai* magmatest gmon.out analysis.txt
