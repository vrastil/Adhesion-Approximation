# Testing
Unit tests are perfomed using [`catch`](https://github.com/catchorg/Catch2).

## Generic tests
For writing generic tests just include at the top of the test file
````
#include <catch.hpp> //< mandatory
#include "test.hpp"  //< optional
````
where ***test.hpp*** contains declarations of test functions used in multiples tests. Tutorial on how to write tests using catch can be found [here](https://github.com/catchorg/Catch2/blob/master/docs/tutorial.md).

## Implementation tests
For testing specific implementation file (`src/*.cpp`) include this file at the top of the test file (to have access to all internal functions) in addition to the above files (***catch.hpp*** and ***test.hpp***). Make sure to name the test file with the same name as the implementation file (including directory) with ***test_*** prefix, e.g. to test file ***src/Foo/Boo.cpp*** create test file ***tests/Foo/test_Boo.cpp***. This is done to avoid multiple definitons errors.