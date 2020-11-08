/*!
 *	\file common.hpp
 *
 *	\author
 *		Andrew Yeliseyev (Russian Technological University, KMBO-03-16, Russia, 2020)
 */
#ifndef _COMMON_HPP_
#define _COMMON_HPP_



#include <cstdint>
#include <cstdio>
#include <functional>
#include <vector>
#include <time.h>





#ifndef TSTESTS_VERBOSITY_LEVEL
#	define TSTESTS_VERBOSITY_LEVEL          0
#endif // TSTESTS_VERBOSITY_LEVEL

#define TSTESTS_TEST_FUNCTION(name)         TsTestsReturnCode const name(TsTestsInfo *const test_info)

#define LOGBLOCK(core)                      do {time(&_time); _ltime = *localtime(&_time); core} while(0);
#define LOGTIMESTAMP                        _ltime.tm_mday, _ltime.tm_mon+1, _ltime.tm_year+1900, _ltime.tm_hour, _ltime.tm_min, _ltime.tm_sec

#if TSTESTS_VERBOSITY_LEVEL > 0
#	define TSTESTS_TEST_FUNCTION_BEGIN(name)        char const *const _test_name = #name; time_t _time; struct tm _ltime; clock_t _start_clock = clock(); TsTestsReturnCode answer = TSTESTS_RETURNCODE_SUCCESS; FILE *_log = (test_info == NULL) ? (stdout) : (test_info->log_file);
#	define TSTESTS_TEST_FUNCTION_END                do {clock_t _finish_clock = clock(); fprintf(_log, "Test completed in %Lf seconds.\n\n", (long double)(_finish_clock - _start_clock) / CLOCKS_PER_SEC); return answer;} while(0);
#else
#	define TSTESTS_TEST_FUNCTION_BEGIN(name)        TsTestsReturnCode answer = TSTESTS_RETURNCODE_SUCCESS;
#	define TSTESTS_TEST_FUNCTION_END                return answer;
#endif // TSTESTS_VERBOSITY_LEVEL 0

#if TSTESTS_VERBOSITY_LEVEL == 1
#	define LOG(message)                             "%s [%02d.%02d.%d %02d:%02d:%02d] %s\n", _test_name, LOGTIMESTAMP, message
#	define PUSHLOG_1(message)                       LOGBLOCK(fprintf(_log, LOG(message));)
#	define LOGF(format, ...)                        "%s [%02d.%02d.%d %02d:%02d:%02d] " format "\n", _test_name, LOGTIMESTAMP, __VA_ARGS__
#	define PUSHLOGF_1(format, ...)                  LOGBLOCK(fprintf(_log, LOGF(format, __VA_ARGS__));)
#else
#	define PUSHLOG_1(message)                       (void) 0;
#	define PUSHLOGF_1(format, ...)                  (void) 0;
#endif // TSTESTS_VERBOSITY_LEVEL 1

#if TSTESTS_VERBOSITY_LEVEL == 2
#	define LOG(message)                             "%s [%02d.%02d.%d %02d:%02d:%02d] %s\n", _test_name, LOGTIMESTAMP, message
#	define PUSHLOG_2(message)                       LOGBLOCK(fprintf(_log, LOG(message));)
#	define LOGF(format, ...)                        "%s [%02d.%02d.%d %02d:%02d:%02d] " format "\n", _test_name, LOGTIMESTAMP, __VA_ARGS__
#	define PUSHLOGF_2(format, ...)                  LOGBLOCK(fprintf(_log, LOGF(format, __VA_ARGS__));)
#else
#	define PUSHLOG_2(message)                       (void) 0;
#	define PUSHLOGF_2(format, ...)                  (void) 0;
#endif // TSTESTS_VERBOSITY_LEVEL 2

#if TSTESTS_VERBOSITY_LEVEL >= 3
#	define LOG(delim, message)                      "%s%s%s\n", _test_name, delim, message
#	define LOGCONT(message)                         "\t%s\n", message
#	define PUSHLOG_3(message)                       LOGBLOCK(fprintf(_log, LOG("\n\t", message));)
#	define APPENDLOG_3(message)                     LOGBLOCK(fprintf(_log, LOGCONT(message));)
#	define LOGF(delim, format, ...)                 "%s%s" format "\n", _test_name, delim, __VA_ARGS__
#	define LOGFCONT(format, ...)                    "\t" format "\n", __VA_ARGS__
#	define PUSHLOGF_3(format, ...)                  LOGBLOCK(fprintf(_log, LOGF("\n\t", format, __VA_ARGS__));)
#	define APPENDLOGF_3(format, ...)                LOGBLOCK(fprintf(_log, LOGFCONT(format, __VA_ARGS__));)
#else
#	define PUSHLOG_3(message)                       (void) 0;
#	define APPENDLOG_3(message)                     (void) 0;
#	define PUSHLOGF_3(format, ...)                  (void) 0;
#	define APPENDLOGF_3(format, ...)                (void) 0;
#endif // TSTESTS_VERBOSITY_LEVEL 3

#if TSTESTS_VERBOSITY_LEVEL >= 4
#	define PUSHLOG_4(message)                       LOGBLOCK(fprintf(_log, LOGF("", " [%02d.%02d.%d %02d:%02d:%02d] %s", LOGTIMESTAMP, message));)
#	define PUSHLOGF_4(format, ...)                  LOGBLOCK(fprintf(_log, LOGF("", " [%02d.%02d.%d %02d:%02d:%02d] " format, LOGTIMESTAMP, __VA_ARGS__));)
#else
#	define PUSHLOG_4(message)                       (void) 0;
#	define PUSHLOGF_4(format, ...)                  (void) 0;
#endif // TSTESTS_VERBOSITY_LEVEL 4



#define TSTESTS_COORDINATE_TYPE     long double
#define TSTESTS_DIGITAL_TYPE        uint64_t

#define TSTESTS_LOG_IN_CONSOLE      stdout



/**
 *  \brief
 *  Stores parameters for \c tstest functions.
 *
 *  Since the set of parameters for all tests is constant, it seems
 *  reasonable to collect all their arguments into one structure and
 *  use it as a sole parameter for test functions.
 *
 *  \param      t                   Assumed defect of the net (to be validated by \c tstest_definition).
 *  \param      m                   Exponent of the amount of points in net (net consists of \f$ 2^m \f$ points).
 *  \param      s                   Dimension.
 *  \param      bitwidth            Number of bits used for points generation (i.e. number of columns in generating matrices).
 *  \param      next_point_getter   A function that returns net's point by its index as a vector of \c TSTEST_DIGITAL_TYPE,
 *                                  if \c TSTESTS_OPTIMISE_FOR_DIGITAL_NETS is defined, or of \c TSTESTS_COORDINATE_TYPE otherwise.
 *  \param      gamma_matrix_getter A function that returns generator's Gamma matrices.
 *  \param      log_file            A valid \c FILE pointer to the opened file or \c TSTESTS_LOG_IN_CONSOLE.
 *
 *  \warning
 *  \b bitwidth parameter is only applicable for digital nets and, thus, will be ignored if
 *  \c TSTESTS_OPTIMISE_FOR_DIGITAL_NETS is not defined.
 *
 *  \warning
 *  \b log_file parameter sets the output destination for logs. It is ignored when
 *  \c TSTESTS_VERBOSITY_LEVEL is set to zero.
 */
typedef struct TsTestsInfo
{
	uint8_t     t;
	uint8_t     m;
	uint8_t     s;
	uint8_t     bitwidth;
#	ifdef TSTESTS_OPTIMISE_FOR_DIGITAL_NETS
	std::function<std::vector<TSTESTS_DIGITAL_TYPE>(uint64_t const)> next_point_getter;
#	else
	std::function<std::vector<TSTESTS_COORDINATE_TYPE>(uint64_t const)> next_point_getter;
#	endif // TSTESTS_OPTIMISE_FOR_DIGITAL_NETS
	std::function<std::vector<std::vector<uint>>(uint64_t const)> gamma_matrix_getter;
	FILE       *log_file;
}
TsTestsInfo;



typedef enum TsTestsReturnCode
{
	TSTESTS_RETURNCODE_SUCCESS              =  0,
	TSTESTS_RETURNCODE_FAIL_GENERAL         = -1,
	TSTESTS_RETURNCODE_FAIL_INPUT           = -2,
	TSTESTS_RETURNCODE_FAIL_MEMORY          = -3
}
TsTestsReturnCode;

#define TSTESTS_PRINT_RETURNCODE(code)      do {                                                                        \
                                                char const *const ret[] = {"TSTESTS_RETURNCODE_FAIL_MEMORY",            \
																		   "TSTESTS_RETURNCODE_FAIL_INPUT",             \
                                                                           "TSTESTS_RETURNCODE_FAIL_GENERAL",           \
                                                                           "TSTESTS_RETURNCODE_SUCCESS"};               \
                                                printf("%s\n", ret[(code) + 3]);                                        \
                                            } while(0);



#endif // _COMMON_HPP_
