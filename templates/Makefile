CXXFLAGS = -std=c++14 -O3 -Wall -Wextra -pedantic -DNDEBUG -I ../../templates/include
LIBS = -ldnest4 -lpthread

default:
	$(CXX) -I ../../../DNest4/code $(CXXFLAGS) -c *.cpp
	$(CXX) -pthread -L ../../../DNest4/code -o main *.o $(LIBS)
	rm *.o