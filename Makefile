montyburns: $(patsubst %.cpp,%.o,$(wildcard *.cpp))
	g++ -ggdb -O3 -lgmp -lrt $^ -o $@

%.o: %.cpp
	g++ -ggdb -O3 -I. -c $< -o $@

clean:
	rm *.o
