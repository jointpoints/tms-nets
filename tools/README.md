### Available `make` options

* `make static_lib` : builds a "tms-nets" static library;
* `make tester` : builds a unit tester for a current version of the library;
* `make docs` : compiles HTML documentation using Doxygen.



### Warning!

Never forget to:

* add all new source files to `UNITS` or `TEST_UNITS` variables in `Makefile` **for both operating systems**;
* change value of `TMS_VERSION` accordingly;
* change value of `TMS_STABILITY` to `stable` before each merge to **master** branch (revert to `unstable` after merge);
* always run `make static_lib` **before** `make tester`, the latter requires a fresh build of the library.



### Tester options

Tester call (`tester.exe` on Windows, `tester` on Linux) may be supplied with additional filters:

* `[<class name>]` : only run tests for the specified class;
* `[nets]` : only run tests for digital nets, equivalent to `[DigitalNet],[Niederreiter],[Sobol]`.

Type the `-s` key after filters (if you enter any) to view the full list of tests.
