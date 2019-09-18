# make
hello : helloworld.o
	@g++ -o hello helloworld.o

helloworld.o : helloworld.cpp
	@g++ -g -c helloworld.cpp

clean :
	@rm *.o hello