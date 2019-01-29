.PHONY: tests

CC=clang++-6.0 #g++

INCLUDE_PATH=${PWD}/include

ifeq ($(DEBUG), 1)
	DEBUG_FLAGS=-g
	OPT_FLAGS=-O0
else
	DEBUG_FLAGS=
	OPT_FLAGS=-O3
endif

ifeq ($(CC), g++)
	OPT_FLAGS+= -fopt-info-vec-optimized
else ifeq ($(CC), clang++-6.0)
	OPT_FLAGS+= -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize
endif

MATH=include/math.h
MATRIX=include/matrix.h
MM_IO=include/mm_io.h utils/mm_io.cpp obj/mm_io.o
UTILS=${MM_IO} include/check_assert.h include/sp_elem_ptr.h

obj/mm_io.o: include/mm_io.h utils/mm_io.cpp
	@mkdir -p obj
	@${CC} ${DEBUG_FLAGS} ${OPT_FLAGS} --std=c++11 -I ${INCLUDE_PATH} -c utils/mm_io.cpp -o obj/mm_io.o

#Tests
tests/test_storage_read_sparse: tests/test_storage_read_sparse.cpp ${MATRIX} ${UTILS}
	@${CC} ${DEBUG_FLAGS} ${OPT_FLAGS} --std=c++11 -I ${INCLUDE_PATH} tests/test_storage_read_sparse.cpp obj/mm_io.o -o tests/test_storage_read_sparse

tests/test_storage_read_dense: tests/test_storage_read_dense.cpp ${MATRIX} ${UTILS}
	@${CC} ${DEBUG_FLAGS} ${OPT_FLAGS} --std=c++11 -I ${INCLUDE_PATH} tests/test_storage_read_dense.cpp obj/mm_io.o -o tests/test_storage_read_dense

tests/test_math_simpleTS: tests/test_math_simpleTS.cpp ${MATH} ${MATRIX} ${UTILS}
	@${CC} ${DEBUG_FLAGS} ${OPT_FLAGS} -fopenmp --std=c++11 -I ${INCLUDE_PATH} tests/test_math_simpleTS.cpp obj/mm_io.o -o tests/test_math_simpleTS

tests: tests/test_storage_read_sparse tests/test_storage_read_dense tests/test_math_simpleTS
