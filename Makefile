
INCLUDE_PATH="${PWD}/include"

ifeq ($(DEBUG), 1)
	DEBUG_FLAGS=-g
else
	DEBUG_FLAGS=
endif

mm_io: include/mm_io.h utils/mm_io.cpp obj/mm_io.o
	g++ ${DEBUG_FLAGS} --std=c++11 -I ${INCLUDE_PATH} -c utils/mm_io.cpp -o obj/mm_io.o

sp_elem_ptr: include/sp_elem_ptr.h

utils: mm_io sp_elem_ptr

include/matrix.h: mm_io

matrix: include/matrix.h


#Tests
tests/test_storage_read_sparse: matrix mm_io
	g++ ${DEBUG_FLAGS} --std=c++11 -I ${INCLUDE_PATH} tests/test_storage_read_sparse.cpp obj/mm_io.o -o tests/test_storage_read_sparse

tests/test_storage_read_dense: matrix mm_io
	g++ ${DEBUG_FLAGS} --std=c++11 -I ${INCLUDE_PATH} tests/test_storage_read_dense.cpp obj/mm_io.o -o tests/test_storage_read_dense

tests: tests/test_storage_read_sparse tests/test_storage_read_dense
