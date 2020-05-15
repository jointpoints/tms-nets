WORKDIR = %cd%

CC = gcc.exe
CXX = g++.exe
AR = ar.exe
LD = g++.exe
WINDRES = windres.exe

INC = 
CFLAGS = -Wall
RESINC = 
LIBDIR = 
LIB = 
LDFLAGS = 

INC_RELEASE = $(INC)
CFLAGS_RELEASE = $(CFLAGS) -O2 -std=c++17
RESINC_RELEASE = $(RESINC)
RCFLAGS_RELEASE = $(RCFLAGS)
LIBDIR_RELEASE = $(LIBDIR)
LIB_RELEASE = $(LIB)
LDFLAGS_RELEASE = $(LDFLAGS) -static
OBJDIR_RELEASE = obj\\Release
DEP_RELEASE = 
OUT_RELEASE = bin\\Release\\tms-nets.exe

OBJ_RELEASE = $(OBJDIR_RELEASE)\\main.o $(OBJDIR_RELEASE)\\tests\\util\\bit_counters.o $(OBJDIR_RELEASE)\\tests\\util\\incremental_pca.o

all: release

clean: clean_release

before_release: 
	cmd /c if not exist bin\\Release md bin\\Release
	cmd /c if not exist $(OBJDIR_RELEASE) md $(OBJDIR_RELEASE)
	cmd /c if not exist $(OBJDIR_RELEASE)\\tests\\util md $(OBJDIR_RELEASE)\\tests\\util

after_release: 

release: before_release out_release after_release

out_release: before_release $(OBJ_RELEASE) $(DEP_RELEASE)
	$(LD) $(LIBDIR_RELEASE) -o $(OUT_RELEASE) $(OBJ_RELEASE)  $(LDFLAGS_RELEASE) $(LIB_RELEASE)

$(OBJDIR_RELEASE)\\main.o: main.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c automatic_tester.cpp -o $(OBJDIR_RELEASE)\\main.o

$(OBJDIR_RELEASE)\\tests\\util\\bit_counters.o: tests\\util\\bit_counters.cpp
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c tests\\util\\bit_counters.cpp -o $(OBJDIR_RELEASE)\\tests\\util\\bit_counters.o

$(OBJDIR_RELEASE)\\tests\\util\\incremental_pca.o: tests\\util\\incremental_pca.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c tests\\util\\incremental_pca.cpp -o $(OBJDIR_RELEASE)\\tests\\util\\incremental_pca.o

clean_release: 
	cmd /c del /f $(OBJ_RELEASE) $(OUT_RELEASE)
	cmd /c rd bin\\Release
	cmd /c rd $(OBJDIR_RELEASE)
	cmd /c rd $(OBJDIR_RELEASE)\\tests\\util

.PHONY: before_release after_release clean_release

