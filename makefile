all: O0 O1 O2 O3 O3M

O0: 
    gcc -o run0 Source.c
	
O1: 
    gcc -O1 -o run1 Source.c
	
O2: 
    gcc -O2 -o run2 Source.c
	
O3: 
    gcc -O3 -o run3 Source.c	

O3M: 
    gcc -O3 -march=native -o run3m Source.c	

clean:
    rm run0 run1 run2 run3 run3m