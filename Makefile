
montyburns: burns.o monty.o main.o primes.o
	g++ -ggdb -O3 burns.o primes.o monty.o main.o -lgmp -lrt -o montyburns

burns.o: burns.h burns.cpp
	g++ -ggdb -O3 -I. -c burns.cpp

primes.o: primes.h primes.cpp
	g++ -ggdb -O3 -I. -c primes.cpp

monty.o: monty.h monty.cpp
	g++ -ggdb -O3 -I. -c monty.cpp

main.o: main.cpp
	g++ -ggdb -O3 -I. -c main.cpp --save-temps

clean:
	rm *.o