.PHONY: tests

INCLUDE_PATH=${PWD}/include

ifeq ($(DEBUG), 1)
	DEBUG_FLAGS=-g
	OPT_FLAGS=-O0
else
	DEBUG_FLAGS=
	OPT_FLAGS=-O3 -fopt-info-vec-optimized
endif

MATH=include/math.h
MATRIX=include/matrix.h
MM_IO=include/mm_io.h utils/mm_io.cpp obj/mm_io.o
UTILS=${MM_IO} include/check_assert.h include/sp_elem_ptr.h

obj/mm_io.o: include/mm_io.h utils/mm_io.cpp
	@mkdir -p obj
	@g++ ${DEBUG_FLAGS} ${OPT_FLAGS} --std=c++11 -I ${INCLUDE_PATH} -c utils/mm_io.cpp -o obj/mm_io.o

#Tests
tests/test_storage_read_sparse: ${MATRIX} ${UTILS}
	@g++ ${DEBUG_FLAGS} ${OPT_FLAGS} --std=c++11 -I ${INCLUDE_PATH} tests/test_storage_read_sparse.cpp obj/mm_io.o -o tests/test_storage_read_sparse

tests/test_storage_read_dense: ${MATRIX} ${UTILS}
	@g++ ${DEBUG_FLAGS} ${OPT_FLAGS} --std=c++11 -I ${INCLUDE_PATH} tests/test_storage_read_dense.cpp obj/mm_io.o -o tests/test_storage_read_dense

tests/test_math_simpleTS: ${MATH} ${MATRIX} ${UTILS}
	@g++ ${DEBUG_FLAGS} ${OPT_FLAGS} -fopenmp --std=c++11 -I ${INCLUDE_PATH} tests/test_math_simpleTS.cpp obj/mm_io.o -o tests/test_math_simpleTS

tests: tests/test_storage_read_sparse tests/test_storage_read_dense tests/test_math_simpleTS
