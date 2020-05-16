CC = gcc.exe
CXX = g++.exe
AR = ar.exe
LD = g++.exe
WINDRES = windres.exe

CFLAGS = -Wall -O2 -std=c++17
LDFLAGS = -static

OBJDIR = automatic_tester_obj
OBJ = $(OBJDIR)\\automatic_tester.o $(OBJDIR)\\util\\bit_counters.o $(OBJDIR)\\util\\incremental_pca.o

OUT = automatic_tester.exe





release: release_prepare release_object1 release_object2 release_object3 release_assemble release_clean

release_prepare:
	powershell New-Item -ItemType Directory -Force -Path $(OBJDIR)\\util

release_object1: automatic_tester.cpp
	$(CXX) $(CFLAGS) $(INC) -c automatic_tester.cpp -o $(OBJDIR)\\automatic_tester.o

release_object2: ..\\util\\bit_counters.cpp
	$(CC) $(CFLAGS) $(INC) -c ..\\util\\bit_counters.cpp -o $(OBJDIR)\\util\\bit_counters.o

release_object3: ..\\util\\incremental_pca.cpp
	$(CXX) $(CFLAGS) $(INC) -c ..\\util\\incremental_pca.cpp -o $(OBJDIR)\\util\\incremental_pca.o

release_assemble:
	$(LD) -o $(OUT) $(OBJ) $(LDFLAGS)

release_clean:
	powershell Remove-Item $(OBJDIR) -Force -Recurse


