
INCLUDE_PATH="${PWD}/include"

mm_io: include/mm_io.h utils/mm_io.cpp
	g++ --std=c++11 -I ${INCLUDE_PATH} -c utils/mm_io.cpp -o obj/mm_io.o

utils: mm_io

lib/matrix.cpp: include/matrix.h
	g++ --std=c++11 -I ${INCLUDE_PATH} -c lib/matrix.cpp -o obj/matrix.o

matrix: include/matrix.h lib/matrix.cpp
