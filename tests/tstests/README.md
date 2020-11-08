TsTests Usage and Development Guide
===================================




## Contents

  * [Purpose of TsTests](#purpose-of-tstests)
  * [Signature](#signature)
  * [Return codes](#return-codes)
  * [Verbosity levels and logs](#verbosity-levels-and-logs)
  * [Optimisations for digital nets](#optimisations-for-digital-nets)
  * [`TsTestsInfo` structure](#tstestsinfo-structure)
  * [Available tests](#available-tests)
    * [1. `tstest_uniqueness`. Test for component-wise uniqueness](#1-tstest_uniqueness-test-for-component-wise-uniqueness)
	  * [1.1. Description](#11-description)
	  * [1.2. Limits of applicability](#12-limits-of-applicability)
	* [2. `tstest_definition`. Test for (t, m, s)-net definition](#2-tstest_definition-test-for-t-m-s-net-definition)
	  * [2.1. Description](#21-description)
	  * [2.2. Limits of applicability](#22-limits-of-applicability)
	* [3. `tstest_principals`. Test for principal components analysis](#3-tstest_principals-test-for-principal-components-analysis)
	  * [3.1. Description](#31-description)
	  * [3.2. Limits of applicability](#32-limits-of-applicability)
  * [Automatic tester](#automatic-tester)
  * [Development guidelines](#development-guidelines)




## Purpose of TsTests

TsTests are created to investigate properties of (t, m, s)-nets. Tests presented here have independent interface from the generator in this repository, so there should be no problems in using them for validating other algorithms. All files that are contained within this folder and its subfolders (except for `README.md` file and `automatic_tester` subfolder) are necessary for the correct operation of tests.

All TsTests can be subdivided into two groups: *validation* tests and *analytical* tests. Validation tests are used to verify particular *hypothesis* for a given test case or, in other words, find the answer for a certain decision problem. Analytical tests, on the other hand, are used to perform some form of calculation and analysis in order to present numeric *characteristics* for the given test case.

[^ to the top ^](#contents)




## Signature

All TsTests have the following signature:

    TsTestsReturnCode const <name of the TsTest>(TsTestsInfo *const test_info)

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

TsTests directly support 5 levels of verbosity that affect the amount of details to be included into logs. Information that is needed to be displayed on each verbosity level varies for validation tests and analytical tests.

  | Verbosity level | Validation tests                                                                                                                                                                                                                                                                                               | Analytical tests                                                                                                                                                             |
  |:---------------:|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
  | `0`             | no logging is performed                                                                                                                                                                                                                                                                                        | no logging is performed                                                                                                                                                      |
  | `1`             | `+`, if `TSTESTS_RETURNCODE_SUCCESS` is returned; `-`, if `TSTESTS_RETURNCODE_FAIL_GENERAL` is returned; `- [!]` otherwise                                                                                                                                                                                     | `+`, if `TSTESTS_RETURNCODE_SUCCESS` is returned; `-`, if `TSTESTS_RETURNCODE_FAIL_GENERAL` is returned; `- [!]` otherwise                                                   |
  | `2`             | `+` and the short statement of verified hypothesis, if `TSTESTS_RETURNCODE_SUCCESS` is returned; `-` and the short statement of unverified hypothesis, if `TSTESTS_RETURNCODE_FAIL_GENERAL` is returned; `- [!]` otherwise                                                                                     | `+`, if `TSTESTS_RETURNCODE_SUCCESS` is returned; `-`, if `TSTESTS_RETURNCODE_FAIL_GENERAL` is returned; `- [!]` otherwise                                                   |
  | `3`             | `Answer: POSITIVE.` and the full statement of verified hypothesis, if `TSTESTS_RETURNCODE_SUCCESS` is returned; `Answer: NEGATIVE.` and the full statement of unverified hypothesis, if `TSTESTS_RETURNCODE_FAIL_GENERAL` is returned; `Answer: NEGATIVE.` and the full description of occured error otherwise | `Answer: POSITIVE.` and all calculated characteristics, if `TSTESTS_RETURNCODE_SUCCESS` is returned; `Answer: NEGATIVE.` and the full description of occured error otherwise |
  | `4`             | logs of all intermediate steps followed by the information described for verbosity level `3`                                                                                                                                                                                                                   | logs of all intermediate steps followed by the information described for verbosity level `3`                                                                                 |

Verbosity level **can** be specified by defining the `TSTESTS_VERBOSITY_LEVEL n` macro where `n` **must** be replaced with one of the values described above. If verbosity level is not manually set, it is implied to be `0`.

[^ to the top ^](#contents)




## Optimisations for digital nets

When testing the digital nets, one **may** apply certain optimisations to make calculations faster and more precise. These optimisations **can** be switched on by defining `TSTESTS_OPTIMISE_FOR_DIGITAL_NETS` macro. Note, that some tests **may** operate for digital nets only, some tests **may** operate for non-digital nets only. These peculiarities are stated in the description of each individual test.

TsTests view every point of a net as a vector of `TSTESTS_COORDINATE_TYPE` floating-point values. However, enabling optimisations changes this to a vector of `TSTESTS_DIGITAL_TYPE` unsigned integral values. Both `TSTESTS_COORDINATE_TYPE` and `TSTESTS_DIGITAL_TYPE` **may** be redefined by user.

[^ to the top ^](#contents)




## `TsTestsInfo` structure

In order to pass parameters into TsTests, the `TsTestsInfo` structure **must** be used. Each structure represents a single test case. `TsTestsInfo` is defined in `util/common.hpp` and contains the following fields:

  * `t` — expected defect of a given (t, m, s)-net (this value is verified by `tstest_definition`, see below for more info);
  * `m` — binary logarithm of the amount of points in a given (t, m, s)-net (net **must** consist of `2^m` points);
  * `s` — dimension of space;
  * `bitwidth` — number of bits used for points generation (meaningful for digital nets only);
  * `next_point_getter` — a function that returns net's point by its index as a vector of `TSTESTS_DIGITAL_TYPE`, if `TSTESTS_OPTIMISE_FOR_DIGITAL_NETS` is defined, or of `TSTESTS_COORDINATE_TYPE` otherwise;
  * `log_file` — `TSTESTS_LOG_IN_CONSOLE` or a valid `FILE` pointer to the opened file.

A pointer to a `TsTestsInfo` structure is the only argument of any TsTest.

[^ to the top ^](#contents)




## Available tests

The following tests are implemented.

#### 1. `tstest_uniqueness`. Test for component-wise uniqueness

###### 1.1. Description

*File*: `tstest_uniqueness.hpp`

*Type*: Validation test.

*Brief*: This test validates component-wise uniqueness of all generated points.

Validation of component-wise uniqueness is performed using bit arrays: one bit array per each dimension. This allows to reduce memory costs and, hence, extend the usage of this test on large-scaled data.

*Returns*:

  * `TSTESTS_RETURNCODE_SUCCESS` in case of successful pass of the test;
  * `TSTESTS_RETURNCODE_FAIL_GENERAL` in case of existence of at least one repetetive coordinate for at least one component;
  * `TSTESTS_RETURNCODE_FAIL_INPUT` in case of invalidity of `test_info` pointer or in case of an attempt to use this test without defined `TSTESTS_OPTIMISE_FOR_DIGITAL_NETS` macro;
  * `TSTESTS_RETURNCODE_FAIL_MEMORY` in case of dynamic memory allocation fail.

###### 1.2. Limits of applicability

| Criterion                         | Value
|-----------------------------------|-------
| Maximum possible `bitwidth` value | 63
| Support of digital nets           | supported
| Support of non-digital nets       | not supported


#### 2. `tstest_definition`. Test for (t, m, s)-net definition

###### 2.1. Description

*File*: `tstest_definition.hpp`

*Type*: Validation test.

*Brief*: This test validates the definition of (t, m, s)-net for the generated set of points.

Validation of definition is performed by calculation of amount of points within each elementary interval. Counters of points occupy the least possible amount of memory due to bitwise packaging.

*Returns*:

  * `TSTESTS_RETURNCODE_SUCCESS` in case of successful definition assertion;
  * `TSTESTS_RETURNCODE_FAIL_GENERAL` in case when generated points fail to meet the requirements of definition by not demonstrating the needed amount of points within at least one elementary interval;
  * `TSTESTS_RETURNCODE_FAIL_INPUT` in case of invalidity of `test_info` pointer or in case when generated points fail to meet the requirements of definition by improperly set `t` and `m` parameters;
  * `TSTESTS_RETURNCODE_FAIL_MEMORY` in case of dynamic memory allocation fail.

###### 2.2. Limits of applicability

| Criterion                         | Value
|-----------------------------------|-------
| Maximum possible `bitwidth` value | 64
| Support of digital nets           | supported
| Support of non-digital nets       | supported


#### 3. `tstest_principals`. Test for principal components analysis

###### 3.1. Description

*File*: `tstest_principals.hpp`

*Type*: Analytical test.

*Brief*: This test performs the incremental principal component analysis for the set of generated points.

This test can be used to find the axes in multidimensional space along which the coordinates of points are variated the most. This is useful for recognising linear patterns in the mutual disposition of points.

*Returns*:

  * `TSTESTS_RETURNCODE_SUCCESS` in case of successful completion of calculations;
  * `TSTESTS_RETURNCODE_FAIL_INPUT` in case of invalidity of `test_info` pointer;
  * `TSTESTS_RETURNCODE_FAIL_MEMORY` in case of dynamic memory allocation fail.

###### 3.2. Limits of applicability

| Criterion                         | Value
|-----------------------------------|-------
| Maximum possible `bitwidth` value | 64
| Support of digital nets           | supported
| Support of non-digital nets       | supported

[^ to the top ^](#contents)




## Automatic tester

`automatic_tester/automatic_tester.cpp` contains an automatic tester for `tms::Niederreiter` generator from this repository. It supports three different modes:

| Mode       | Description                                                                                 | Approximate runtime |
|:----------:|---------------------------------------------------------------------------------------------|:-------------------:|
| Critical   | A set of 5 simple test cases covering the most obvious potential vulnerabilities of program | 1~5 seconds         |
| Regular    | A set of 45 test cases for various (t, m, s)-nets with *s* varying from 1 to 9              | 1~5 minutes         |
| Exhaustive | A single test case to validate the results for large-scale nets                             | 3~5 hours           |

Tester **can** be built with the help of `automatic_tester/automatic_tester.mak` using the following command line:

    make -f automatic_tester.mak

Tester **can** be supplied with the following command line arguments:

  * `-c`, `-r` or `-e` — specifies the mode (critical, regular or exhaustive, respectively); critical tests will be launched when these keys are omitted.
  * `-log f` — specifies file to write logs by its path `f`; log will be performed into console when this argument is omitted.

You **may** look through `automatic_tester/automatic_tester_log_*.txt` files to see the results for the current version of generator.

[^ to the top ^](#contents)




## Development guidelines

Any developer who wishes to design a new TsTest **must** follow these guidelines in order to organically fit into their infrastructure.

  * Any TsTest **must** be defined by `TSTESTS_TEST_FUNCTION(name)` where `name` **must** be replaced with the desired name of test function.
  * The first line of code in the body of any TsTest **must** be `TSTESTS_TEST_FUNCTION_BEGIN(name)` where `name` **must** be replaced with capitalised name of test function.
  * Name of any TsTest **should** begin with `tstest_`.
  * A pointer to a `TsTestsInfo` structure is the only argument of any TsTest. This pointer **can** be accessed in the body of any TsTest under the name of `test_info`.
  * Any TsTest **should** be explicitly intolerant towards invalid `test_info` pointer.
  * Each TsTest is automatically supplied with the variable `answer` of a `TsTestsReturnCode` type. Its default value is `TSTESTS_RETURNCODE_SUCCESS` which **can** be freely changed in the main body of TsTest.
  * Any TsTest **must** use `TSTESTS_TEST_FUNCTION_END` instead of `return ...;` statements.
  * `TSTESTS_TEST_FUNCTION_END` finishes execution of any TsTest and returns the value of the `answer` variable. Thus, the final result of any TsTest **must** be stored within the `answer` variable before the execution of any of `TSTESTS_TEST_FUNCTION_END`.
  * Return values of any TsTest **should** be construed in accordance with the [Return codes](#return-codes) section of this document.
  * Any TsTest **must** provide logging for verbosity levels mentioned in the [Verbosity levels and logs](#verbosity-levels-and-logs) section of this document.
  * If any TsTest returns value in complete accordance with the [Return codes](#return-codes) section of this document, then its logging **must** be done in accordance with the [Verbosity levels and logs](#verbosity-levels-and-logs) section of this document.
  * Any TsTest **should** tend to use native logging abilities of TsTests infrastructure, specifically `PUSHLOGn(message)`, `PUSHLOGFn(format, ...)`, `APPENDLOG3(message)` and `APPENDLOGF3(format, ...)` macros where `n` **must** be replaced with corresponding verbosity level from `1` to `4`.
  * Any TsTest **should** be able to work with both digital and non-digital nets.
  * Any TsTest **must** be contained within a single `.hpp` file separately from other TsTests and any other unrelated code.

[^ to the top ^](#contents)