## Building the tests

To build the driver binaries,
```sh
$ make tests
```

To test the sparse-triangular-solve routines,
```sh
$ tests/test_math_simpleTS use_double path/to/sparse/A path/to/dense/B
```
|use_double|type of A[i][j]|
|----------|:--------------|
|0         |`float`        |
|1         |`double`       |


