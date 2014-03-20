test.o: community.o data.o
	g++ community.o data.o main.cpp -o lotto

data.o: data.cpp data.h community.o
	g++ -c data.cpp -o data.o

community.o : community.cpp community.h 
	g++ -c community.cpp -o community.o

clean:
	rm data.o community.o
