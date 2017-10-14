SHELL:=/bin/bash
all: defgcc p11 p12 p13 p23 p24O0 p24O1 p24O2 p24O3 newgcc p24N0 p24N1 p24N2 p24N3
	
defgcc:
	module purge

p11:
	gcc -o p11 1p1.c

p12:
	gcc -o p12 1p2.c

p13:
	gcc -o p13 1p3.c

p23:
	gcc -lrt -o p23 2p3.c
	
p24O0:
	gcc -lrt -o p24O0 2p4.c
	
p24O1:
	gcc -lrt -O1 -o p24O1 2p4.c

p24O2:
	gcc -lrt -O2 -o p24O2 2p4.c

p24O3:
	gcc -lrt -O3 -o p24O3 2p4.c

newgcc:
	module load gcc-4.7.2

p24N0:
	gcc -lrt -o p24N0 2p4.c
	
p24N1:
	gcc -lrt -O1 -o p24N1 2p4.c

p24N2:
	gcc -lrt -O2 -o p24N2 2p4.c

p24N3:
	gcc -lrt -O3 -o p24N3 2p4.c

clean:
	rm p11 p12 p13 p23 p24O0 p24O1 p24O2 p24O3 p24N0 p24N1 p24N2 p24N3
	
run:
	qsub p11.job
	qsub p12.job
	qsub p13.job
	qsub p23.job
	qsub p24O0.job
	qsub p24O1.job
	qsub p24O2.job
	qsub p24O3.job
	qsub p24N0.job
	qsub p24N1.job
	qsub p24N2.job
	qsub p24N3.job
	