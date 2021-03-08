### Available `make` options

* `make static_lib` : builds a "tms-nets" static library;
* `make docs` : compiles HTML documentation using Doxygen.



### Warning!

Never forget to:

* add all new source files to `UNITS` variable in `Makefile` **for both operating systems**;
* change value of `TMS_VERSION` accordingly;
* change value of `TMS_STABILITY` to `stable` before each merge to **master** branch (revert to `unstable` after merge).
