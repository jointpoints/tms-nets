TsTests Usage and Development Guide
===================================




## Contents

  * [Purpose of TsTests](#purpose-of-tstests)
  * [Return codes](#return-codes)
  * [Verbosity levels and logs](#verbosity-levels-and-logs)
  * [Optimisations for digital nets](#optimisations-for-digital-nets)
  * [`TsTestsInfo` structure](#tstestsinfo-structure)
  * [Available tests](#available-tests)
    * [1. `tstest_uniqueness`. Test for component-wise uniqueness](#1-tstest-uniqueness-test-for-component-wise-uniqueness`)
	  * [1.1. Description](#1-1-description)
	  * [1.2. Limits of applicability](#1-2-limits-of-applicability)
	* [2. `tstest_definition`. Test for (t, m, s)-net definition](#2-tstest-definition-test-for-t-m-s-net-definition)
	  * [2.1. Description](#2-1-description)
	  * [2.2. Limits of applicability](#2-2-limits-of-applicability)
	* [3. `tstest_principals`. Test for principal components analysis](#3-tstest-principals-test-for-principal-components-analysis)
	  * [3.1. Description](#3-1-description)
	  * [3.2. Limits of applicability](#3-2-limits-of-applicability)
  * [Automatic tester](#automatic-tester)
  * [Development guidelines](#development-guidelines)




## Purpose of TsTests

TsTests are created to verify properties of (t, m, s)-nets constructed with the help of (t, s)-sequences. Tests presented here have independent interface from the generator in this repository, so there should be no problems in using them for validating other programs. All files that are contained within this folder and its subfolders (except for `automatic_tester.mak`, `automatic_tester.cpp`, `automatic_tester_log.txt` and `README.md`) are necessary for the correct operation of tests.
All TsTests can be subdivided into two groups: *validation* tests and *analytical* tests. Validation tests are used to verify particular *hypothesis* for a given test case or, in other words, find the answer for a certain yes/no question. Analytical tests, on the other hand, are used to perform some form of calculation and analysis in order to present numeric *characteristics* for the given test case.
TsTests view every point of net as a vector of `TSTESTS_COORDINATE_TYPE` values. `TSTESTS_COORDINATE_TYPE` **may** be redefined by user.

[^ to the top ^](#contents)




## Return codes

TsTests use the following return codes defined in `util/common.hpp::TsTestsReturnCode`.

  | Return code                       | Validation tests                                                                                       | Analytical tests                                                                               |
  |:---------------------------------:|--------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------|
  | `TSTESTS_RETURNCODE_SUCCESS`      | supposed to be used when test finishes its work successfully and tested hypothesis is found to be true | supposed to be used when test finishes its work successfully and all characteristics are found |
  | `TSTESTS_RETURNCODE_FAIL_GENERAL` | supposed to be used when a contradiction with the tested hypothesis has been discovered                | supposed to be used when test is not able to calculate at least one characteristic             |
  | `TSTESTS_RETURNCODE_FAIL_INPUT`   | supposed to be used when input data is invalid                                                         | supposed to be used when input data is invalid                                                 |
  | `TSTESTS_RETURNCODE_FAIL_MEMORY`  | supposed to be used when memory allocation fails                                                       | supposed to be used when memory allocation fails                                               |

Each return code `r` **may** be printed into console in a form of human-readable string using `TSTESTS_PRINT_RETURNCODE(r)`.

[^ to the top ^](#contents)




## Verbosity levels and logs

TsTests directly support 5 levels of verbosity that affect the amount of details to be included into logs.

  * `0` — no logging is done;
  * `1` — short answers and work time are supposed to be logged;
  * `2` — more detailed answers and work time are supposed to be logged;
  * `3` — complete answers and work time are supposed to be logged;
  * `4` — intermediate steps, complete answers and work time are supposed to be logged.

Verbosity level **can** be determined by defining the `TSTESTS_VERBOSITY_LEVEL n` macro where `n` **must** be replaced with one of the values described above. If verbosity level is not manually set, it is implied to be `0`.

[^ to the top ^](#contents)




## Optimisations for digital nets

When testing the digital nets, one **may** apply certain optimisations to make calculations faster and more precise. These optimisations **can** be switched on by defining `TSTESTS_OPTIMISE_FOR_DIGITAL_NETS` macro. Note, that some tests **may** operate for digital nets only, some tests **may** operate for non-digital nets only. These pacularities are stated in the description of each individual test.

[^ to the top ^](#contents)




## `TsTestsInfo` structure

In order to pass parameters into TsTests, the `TsTestsInfo` structure **must** be used. Each structure represents a single test case. `TsTestsInfo` is defined in `util/common.hpp` and contains the following fields:

  * `t` — expected defect of a given (t, m, s)-net (this value is verified by `tstest_definition`, see below for more info);
  * `m` — binary logarithm of the amount of points in a given (t, m, s)-net (net **must** consist of `2^m` points);
  * `s` — dimension of space;
  * `bitwidth` — number of bits used for points generation (applicable for digital nets only);
  * `next_point_getter` — a function that returns net's point by its index as a vector of `TSTESTS_COORDINATE_TYPE`;
  * `log_file` — `TSTESTS_LOG_IN_CONSOLE` or a valid `FILE` pointer to the opened file.

[^ to the top ^](#contents)




## Available tests

The following tests are implemented.

#### 1. `tstest_uniqueness`. Test for component-wise uniqueness

###### 1.1. Description

*File*: `tstest_uniqueness.hpp`

*Signature*: `TsTestsReturnCode const tstest_uniqueness(TsTestsInfo *const test_info)`

*Type*: Validation test.

*Brief*: This test validates component-wise uniqueness of all generated points.

Validation of component-wise uniqueness is performed using bit arrays: one bit array per each dimension. This allows reducing memory costs and, hence, extends the usage of this test on large-scaled data.

*Returns*:

  * `TSTESTS_RETURNCODE_SUCCESS` in case of successful pass of the test;
  * `TSTESTS_RETURNCODE_FAIL_GENERAL` in case of existence of at least one repetetive coordinate for at least one component;
  * `TSTESTS_RETURNCODE_FAIL_INPUT` in case of invalidity of `test_info` pointer;
  * `TSTESTS_RETURNCODE_FAIL_MEMORY` in case of dynamic memory allocation fail.

###### 1.2. Limits of applicability

| Criterion                        | Value
|----------------------------------|-------
| Maximum advised `bitwidth` value | 32
| Support of digital nets          | supported
| Support of non-digital nets      | not supported


#### 2. `tstest_definition`. Test for (t, m, s)-net definition

###### 2.1. Description

*File*: `tstest_definition.hpp`

*Signature*: `TsTestsReturnCode const tstest_definition(TsTestsInfo *const test_info)`

*Type*: Validation test.

*Brief*: This test validates the definition of (t, m, s)-net for the generated set of points.

Validation of definition is performed by calculation of points within each elementary interval. Counters of points occupy the least possible amount of memory using the bitwise packaging.

*Returns*:

  * `TSTESTS_RETURNCODE_SUCCESS` in case of successful definition assertion;
  * `TSTESTS_RETURNCODE_FAIL_GENERAL` in case when generated points fail to meet the requirements of definition by not demonstrating the needed amount of points within at least one elementary interval;
  * `TSTESTS_RETURNCODE_FAIL_INPUT` in case of invalidity of `test_info` pointer or in case when generated points fail to meet the requirements of definition by improperly set `t` and `m` parameters;
  * `TSTESTS_RETURNCODE_FAIL_MEMORY` in case of dynamic memory allocation fail.

###### 2.2. Limits of applicability

| Criterion                        | Value
|----------------------------------|-------
| Maximum advised `bitwidth` value | 63
| Support of digital nets          | supported
| Support of non-digital nets      | supported


#### 3. `tstest_principals`. Test for principal components analysis

###### 3.1. Description

*File*: `tstest_principals.hpp`

*Signature*: `TsTestsReturnCode const tstest_principals(TsTestsInfo *const test_info)`

*Type*: Analytical test.

*Brief*: This test performs the incremental principal component analysis for the set of generated points.

This test can be used to find the axes in multidimensional space along which the coordinates of points are variated the most. This is useful for recognising linear patterns in the mutual disposition of points.

*Returns*:

  * `TSTESTS_RETURNCODE_SUCCESS` in case of successful completion of calculations;
  * `TSTESTS_RETURNCODE_FAIL_INPUT` in case of invalidity of `test_info` pointer;
  * `TSTESTS_RETURNCODE_FAIL_MEMORY` in case of dynamic memory allocation fail.

###### 3.2. Limits of applicability

| Criterion                        | Value
|----------------------------------|-------
| Maximum advised `bitwidth` value | 64
| Support of digital nets          | supported
| Support of non-digital nets      | supported

[^ to the top ^](#contents)




## Automatic tester

`automatic_tester.cpp` contains an automatic tester for `sequences::Niederreiter` generator from this repository. This tester **can** be built with the help of `automatic_tester.mak` Make-file, if one uses GCC compilers. Tester **can** be supplied with the following command line arguments:

  * `-log f` — specifies file to write logs by its path `f`; log will be performed into console when this argument is omitted.

You **may** look through `automatic_tester_log.txt` to see its output.

[^ to the top ^](#contents)




## Development guidelines

Any developer who wishes to design a new TsTest **must** follow these guidelines in order to organically fit into their infrastructure.

  * Any TsTest **must** return `TsTestsReturnCode const` value in accordance with the [Return codes](#return-codes) section of this document.
  * Any TsTest **must** accept `TsTestsInfo *const` as its only argument.
  * Name of any TsTest **should** begin with `tstest_`.
  * The first line of code in any TsTest **must** be `TSTESTS_TEST_FUNCTION_BEGIN(name, log_file)` where `name` **must** be replaced with capitalised name of test function and `log_file` **must** be replaced with `TSTESTS_LOG_IN_CONSOLE` or a valid pointer to an opened `FILE`.
  * The line of code in any TsTest before any `return` **must** be `TSTESTS_TEST_FUNCTION_END`.
  * Any TsTest **should** support logging and verbosity in accordance with the [Verbosity levels and logs](#verbosity-levels-and-logs) section of this document.
  * Any TsTest **should** tend to use native logging abilities of TsTests infrastructure, specifically `PUSHLOGn(message)`, `PUSHLOGFn(format, ...)`, `APPENDLOG3(message)` and `APPENDLOGF3(format, ...)` macros where `n` **must** be replaced with corresponding verbosity level from `1` to `4`.
  * Any TsTest **should** be able to work with both digital and non-digital nets.

[^ to the top ^](#contents)