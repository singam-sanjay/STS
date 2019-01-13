
INCLUDE_PATH="${PWD}/include"

lib/matrix.cpp: include/matrix.h
	g++ --std=c++11 -I ${INCLUDE_PATH} -c lib/matrix.cpp -o obj/matrix.o

matrix: include/matrix.h lib/matrix.cpp
