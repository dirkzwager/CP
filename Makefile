ca: main.o CA.o Cell.o Wire.o Genome.o Vector.o
	g++ -o ca main.o CA.o Cell.o Wire.o Genome.o Vector.o -lreadline -lglut -lpthread

main.o: main.cpp
	g++ -c main.cpp

CA.o: CA.cpp CA.h
	g++ -c CA.cpp

Cell.o: Cell.cpp Cell.h
	g++ -c Cell.cpp

Wire.o: Wire.cpp Wire.h
	g++ -c Wire.cpp

Genome.o: Genome.cpp Genome.h
	g++ -c Genome.cpp

Vector.o: Vector.cpp Vector.h
	g++ -c Vector.cpp

clean:
	rm -rf *.o