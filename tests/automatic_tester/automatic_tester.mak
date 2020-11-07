CFLAGS = -Wall -O2 -std=c++17
OBJDIR = automatic_tester_obj


ifeq ($(OS),Windows_NT)

CC = gcc.exe
CXX = g++.exe
AR = ar.exe
LD = g++.exe
LDFLAGS = -static
WINDRES = windres.exe
OBJ = $(OBJDIR)\\automatic_tester.o $(OBJDIR)\\util\\bit_counters.o $(OBJDIR)\\util\\incremental_pca.o $(OBJDIR)\\util\\raref.o
OUT = automatic_tester.exe
	
release: release_prepare_win release_object1_win release_object2_win release_object3_win release_object4_win release_assemble_win release_clean_win

else

ifeq ($(shell uname -s), Linux)
LDFLAGS = -static -pthread
else
LDFLAGS = -pthread
endif

CC = gcc
CXX = g++
LD = g++
OBJ = $(OBJDIR)/automatic_tester.o $(OBJDIR)/util/bit_counters.o $(OBJDIR)/util/incremental_pca.o $(OBJDIR)/util/raref.o 
OUT = automatic_tester.out

release: release_prepare_unix release_object1_unix release_object2_unix release_object3_unix release_object4_unix release_assemble_unix release_clean_unix

endif

#TARGETS FOR WINDOWS

release_prepare_win:
	powershell New-Item -ItemType Directory -Force -Path $(OBJDIR)\\util

release_object1_win: automatic_tester.cpp
	$(CXX) $(CFLAGS) $(INC) -c automatic_tester.cpp -o $(OBJDIR)\\automatic_tester.o

release_object2_win: ..\\util\\bit_counters.cpp
	$(CC) $(CFLAGS) $(INC) -c ..\\util\\bit_counters.cpp -o $(OBJDIR)\\util\\bit_counters.o

release_object3_win: ..\\util\\incremental_pca.cpp
	$(CXX) $(CFLAGS) $(INC) -c ..\\util\\incremental_pca.cpp -o $(OBJDIR)\\util\\incremental_pca.o

release_object4_win: ..\\util\\raref.cpp
	$(CXX) $(CFLAGS) $(INC) -c ..\\util\\raref.cpp -o $(OBJDIR)\\util\\raref_pca.o

release_assemble_win:
	$(LD) -o $(OUT) $(OBJ) $(LDFLAGS)

release_clean_win:
	powershell Remove-Item $(OBJDIR) -Force -Recurse


#TARGETS FOR LINUX/OSX

release_prepare_unix:
	if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi; \
	if [ ! -d $(OBJDIR)/util ]; then mkdir $(OBJDIR)/util; fi;

release_object1_unix: automatic_tester.cpp
	$(CXX) $(CFLAGS) $(INC) -c automatic_tester.cpp -o $(OBJDIR)/automatic_tester.o

release_object2_unix: ../util/bit_counters.cpp
	$(CC) $(CFLAGS) $(INC) -c ../util/bit_counters.cpp -o $(OBJDIR)/util/bit_counters.o

release_object3_unix: ../util/incremental_pca.cpp
	$(CXX) $(CFLAGS) $(INC) -c ../util/incremental_pca.cpp -o $(OBJDIR)/util/incremental_pca.o

release_object4_unix: ../util/raref.cpp
	$(CXX) $(CFLAGS) $(INC) -c ../util/raref.cpp -o $(OBJDIR)/util/raref.o

release_assemble_unix:
	$(LD) -o $(OUT) $(OBJ) $(LDFLAGS)

release_clean_unix:
	rm -rf $(OBJDIR)
