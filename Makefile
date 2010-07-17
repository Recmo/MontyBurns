FLAGS = -ggdb -O3 -I.

montyburns: montyburns.o burns.o monty.o main.o primes.o
	g++ ${FLAGS} montyburns.o burns.o primes.o monty.o main.o -lgmp -lrt -o montyburns

montyburns.o: primes.h monty.h burns.h montyburns.h montyburns.cpp
	g++ ${FLAGS} -c montyburns.cpp

burns.o: primes.h monty.h burns.h burns.cpp
	g++ ${FLAGS} -c burns.cpp

primes.o: primes.h primes.cpp
	g++ ${FLAGS} -c primes.cpp

monty.o: primes.h monty.h monty.cpp
	g++ ${FLAGS} -c monty.cpp

main.o: burns.h monty.h primes.h main.cpp
	g++ ${FLAGS} -c main.cpp

clean:
	rm *.o
