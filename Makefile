
INCLUDE_PATH="${PWD}/include"

ifeq ($(DEBUG), 1)
	DEBUG_FLAGS=-g
	OPT_FLAGS=-O0
else
	DEBUG_FLAGS=
	OPT_FLAGS=-O3
endif

mm_io: include/mm_io.h utils/mm_io.cpp obj/mm_io.o
	g++ ${DEBUG_FLAGS} ${OPT_FLAGS} --std=c++11 -I ${INCLUDE_PATH} -c utils/mm_io.cpp -o obj/mm_io.o

include/sp_elem_ptr.h:

include/sp_elem_ptr.h:

include/check_utils.h:

utils: mm_io include/sp_elem_ptr.h include/check_utils.h

include/matrix.h: mm_io

matrix: include/matrix.h

include/math.h: matrix utils

math: include/math.h


#Tests
tests/test_storage_read_sparse: matrix utils
	g++ ${DEBUG_FLAGS} ${OPT_FLAGS} --std=c++11 -I ${INCLUDE_PATH} tests/test_storage_read_sparse.cpp obj/mm_io.o -o tests/test_storage_read_sparse

tests/test_storage_read_dense: matrix utils
	g++ ${DEBUG_FLAGS} ${OPT_FLAGS} --std=c++11 -I ${INCLUDE_PATH} tests/test_storage_read_dense.cpp obj/mm_io.o -o tests/test_storage_read_dense

tests/test_math_simpleTS: matrix math utils
	g++ ${DEBUG_FLAGS} ${OPT_FLAGS} -fopt-info-vec-optimized -fopenmp --std=c++11 -I ${INCLUDE_PATH} tests/test_math_simpleTS.cpp obj/mm_io.o -o tests/test_math_simpleTS

tests: tests/test_storage_read_sparse tests/test_storage_read_dense tests/test_math_simpleTS
