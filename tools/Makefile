TMS_VERSION = 3.0.0
# Change TMS_STABILITY to "stable" before each merge with the "master" branch!
TMS_STABILITY = stable

CPP_COMPILER = g++
ARCHIVER = ar

COMPILER_FLAGS = -O2 -std=c++17
OBJECT_FOLDER = tms_nets_obj

STATIC_LIB_FOLDER = tms-nets (static library, v.$(TMS_VERSION), $(TMS_STABILITY))
STATIC_LIB_OUT_FILE = tms-nets.a





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

static_lib: static_prepare_win $(UNITS) static_assemble_win static_clean_win

docs:
	powershell Set-Location ..; doxygen .\\tools\\doxygen_config; powershell Set-Location .\\tools

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

static_lib: static_prepare_unix $(UNITS) static_assemble_unix static_clean_unix

docs:
	cd ..; doxygen ./tools/doxygen_config; cd ./tools

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
