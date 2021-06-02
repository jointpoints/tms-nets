TMS_VERSION = 3.0.1
# Change TMS_STABILITY to "stable" before each merge with the "master" branch!
TMS_STABILITY = unstable

CPP_COMPILER = g++
ARCHIVER = ar
LINKER = g++

COMPILER_FLAGS = -O2 -std=c++17
OBJECT_FOLDER = tms_nets_obj

STATIC_LIB_FOLDER = tms-nets (static library, v.$(TMS_VERSION), $(TMS_STABILITY))
STATIC_LIB_OUT_FILE = tms-nets.a

TESTER_OBJECT_FOLDER = tester_obj
TESTER_OUT_FILE = tester





ifeq ($(OS),Windows_NT)

# "..\\" before "source" is omitted in SOURCE_FOLDER due to erroneous interpretation by make; it is manually added where needed
SOURCE_FOLDER = source
UNITS = $(SOURCE_FOLDER)\\thirdparty\\irrpoly\\gf.cpp $(SOURCE_FOLDER)\\thirdparty\\irrpoly\\gfpoly.cpp $(SOURCE_FOLDER)\\thirdparty\\irrpoly\\gfcheck.cpp\
        $(SOURCE_FOLDER)\\details\\common.cpp $(SOURCE_FOLDER)\\details\\gf2poly.cpp $(SOURCE_FOLDER)\\details\\genmat.cpp $(SOURCE_FOLDER)\\details\\recseq.cpp\
        $(SOURCE_FOLDER)\\digital_net.cpp $(SOURCE_FOLDER)\\niederreiter.cpp $(SOURCE_FOLDER)\\modified_niederreiter.cpp\
        $(SOURCE_FOLDER)\\analysis\\t.cpp $(SOURCE_FOLDER)\\analysis\\scatter_defect.cpp

INCLUDE_FOLDER = ..\\include
LICENSE_TMS_FILE = ..\\LICENSE.md
LICENSE_TMS_FILE_OUT = LICENSE.txt

TEST_FOLDER = tests
TEST_UNITS_FOLDER = $(TEST_FOLDER)\\units
TEST_UNITS = $(TEST_FOLDER)\\catch2\\catch_amalgamated.cpp $(TEST_FOLDER)\\unit_tests.cpp\
             $(TEST_UNITS_FOLDER)\\unit_DigitalNet.cpp $(TEST_UNITS_FOLDER)\\unit_Niederreiter.cpp

static_lib: static_prepare_win $(UNITS) static_assemble_win static_clean_win

docs:
	powershell Set-Location ..; doxygen .\\tools\\doxygen_config; powershell Set-Location .\\tools

tester: tester_prepare_win $(TEST_UNITS) tester_assemble_win tester_clean_win

else

# "../" before "source" is omitted in SOURCE_FOLDER due to erroneous interpretation by make; it is manually added where needed
SOURCE_FOLDER = source
UNITS = $(SOURCE_FOLDER)/thirdparty/irrpoly/gf.cpp $(SOURCE_FOLDER)/thirdparty/irrpoly/gfpoly.cpp $(SOURCE_FOLDER)/thirdparty/irrpoly/gfcheck.cpp\
        $(SOURCE_FOLDER)/details/common.cpp $(SOURCE_FOLDER)/details/gf2poly.cpp $(SOURCE_FOLDER)/details/genmat.cpp $(SOURCE_FOLDER)/details/recseq.cpp\
        $(SOURCE_FOLDER)/digital_net.cpp $(SOURCE_FOLDER)/niederreiter.cpp $(SOURCE_FOLDER)/modified_niederreiter.cpp\
        $(SOURCE_FOLDER)/analysis/t.cpp $(SOURCE_FOLDER)/analysis/scatter_defect.cpp

INCLUDE_FOLDER = ../include
LICENSE_TMS_FILE = ../LICENSE.md
LICENSE_TMS_FILE_OUT = LICENSE.txt

TEST_FOLDER = tests
TEST_UNITS_FOLDER = $(TEST_FOLDER)/units
TEST_UNITS = $(TEST_FOLDER)/catch2/catch_amalgamated.cpp $(TEST_FOLDER)/unit_tests.cpp\
             $(TEST_UNITS_FOLDER)/unit_DigitalNet.cpp $(TEST_UNITS_FOLDER)/unit_Niederreiter.cpp

static_lib: static_prepare_unix $(UNITS) static_assemble_unix static_clean_unix

docs:
	cd ..; doxygen ./tools/doxygen_config; cd ./tools

tester: tester_prepare_unix $(TEST_UNITS) tester_assemble_unix tester_clean_unix

endif





# STATIC_LIB TARGETS FOR WINDOWS

static_prepare_win:
	powershell New-Item -ItemType Directory -Force -Path $(OBJECT_FOLDER)
	powershell New-Item -ItemType Directory -Force -Path \"$(STATIC_LIB_FOLDER)\"
	powershell Copy-Item -Path $(INCLUDE_FOLDER)\\tms-nets -Destination \"$(STATIC_LIB_FOLDER)\\tms-nets\" -Recurse
	powershell Remove-Item \"$(STATIC_LIB_FOLDER)\\tms-nets\\thirdparty\\jacobi\" -Force -Recurse
	powershell Copy-Item -Path $(INCLUDE_FOLDER)\\tms-nets.hpp -Destination \"$(STATIC_LIB_FOLDER)\\tms-nets.hpp\"
	powershell Copy-Item -Path $(LICENSE_TMS_FILE) -Destination \"$(STATIC_LIB_FOLDER)\\$(LICENSE_TMS_FILE_OUT)\"

$(SOURCE_FOLDER)\\%.cpp:
	$(CPP_COMPILER) $(COMPILER_FLAGS) -c $(addprefix ..\\,$@) -o $(OBJECT_FOLDER)\\$(addsuffix .o,$(basename $(notdir $@)))

static_assemble_win:
	$(ARCHIVER) rvs "$(STATIC_LIB_FOLDER)\\$(STATIC_LIB_OUT_FILE)" $(addprefix $(OBJECT_FOLDER)\\,$(addsuffix .o,$(basename $(notdir $(UNITS)))))

static_clean_win:
	powershell Remove-Item $(OBJECT_FOLDER) -Force -Recurse



# STATIC_LIB TARGETS FOR LINUX/OSX

static_prepare_unix:
	mkdir -p $(OBJECT_FOLDER)
	mkdir -p "$(STATIC_LIB_FOLDER)"
	cp -r $(INCLUDE_FOLDER)/tms-nets "$(STATIC_LIB_FOLDER)/tms-nets"
	rm -rf "$(STATIC_LIB_FOLDER)/tms-nets/thirdparty/jacobi"
	cp $(INCLUDE_FOLDER)/tms-nets.hpp "$(STATIC_LIB_FOLDER)/tms-nets.hpp"
	cp $(LICENSE_TMS_FILE) "$(STATIC_LIB_FOLDER)/$(LICENSE_TMS_FILE_OUT)"

$(SOURCE_FOLDER)/%.cpp:
	$(CPP_COMPILER) $(COMPILER_FLAGS) -c $(addprefix ../,$@) -o $(OBJECT_FOLDER)/$(addsuffix .o,$(basename $(notdir $@)))

static_assemble_unix:
	$(ARCHIVER) rvs "$(STATIC_LIB_FOLDER)/$(STATIC_LIB_OUT_FILE)" $(addprefix $(OBJECT_FOLDER)/,$(addsuffix .o,$(basename $(notdir $(UNITS)))))

static_clean_unix:
	rm -rf $(OBJECT_FOLDER)





# TESTER TARGETS FOR WINDOWS

tester_prepare_win:
	powershell New-Item -ItemType Directory -Force -Path $(TESTER_OBJECT_FOLDER)

$(TEST_UNITS_FOLDER)\\%.cpp:
	$(CPP_COMPILER) $(COMPILER_FLAGS) -c $(addprefix ..\\,$@) -o $(TESTER_OBJECT_FOLDER)\\$(addsuffix .o,$(basename $(notdir $@)))

$(TEST_FOLDER)\\unit_tests.cpp $(TEST_FOLDER)\\catch2\\catch_amalgamated.cpp:
	$(CPP_COMPILER) $(COMPILER_FLAGS) -c $(addprefix ..\\,$@) -o $(TESTER_OBJECT_FOLDER)\\$(addsuffix .o,$(basename $(notdir $@)))

tester_assemble_win:
	$(LINKER) -o $(TESTER_OUT_FILE).exe $(addprefix $(TESTER_OBJECT_FOLDER)\\,$(addsuffix .o,$(basename $(notdir $(TEST_UNITS))))) "$(STATIC_LIB_FOLDER)\\$(STATIC_LIB_OUT_FILE)"

tester_clean_win:
	powershell Remove-Item $(TESTER_OBJECT_FOLDER) -Force -Recurse



# TESTER TARGETS FOR LINUX/OSX

tester_prepare_unix:
	mkdir -p $(TESTER_OBJECT_FOLDER)

$(TEST_UNITS_FOLDER)/%.cpp:
	$(CPP_COMPILER) $(COMPILER_FLAGS) -c $(addprefix ../,$@) -o $(TESTER_OBJECT_FOLDER)/$(addsuffix .o,$(basename $(notdir $@)))

$(TEST_FOLDER)/unit_tests.cpp $(TEST_FOLDER)/catch2/catch_amalgamated.cpp:
	$(CPP_COMPILER) $(COMPILER_FLAGS) -c $(addprefix ../,$@) -o $(TESTER_OBJECT_FOLDER)/$(addsuffix .o,$(basename $(notdir $@)))

tester_assemble_unix:
	$(LINKER) -o $(TESTER_OUT_FILE) $(addprefix $(TESTER_OBJECT_FOLDER)/,$(addsuffix .o,$(basename $(notdir $(TEST_UNITS))))) "$(STATIC_LIB_FOLDER)/$(STATIC_LIB_OUT_FILE)"

tester_clean_unix:
	rm -rf  $(TESTER_OBJECT_FOLDER)
