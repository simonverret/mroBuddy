OPTIONS = -O2 -Wall
EXEC = mro
INSTALLDIR := $(HOME)/bin/
#COMPILER = gcc
COMPILER = /Users/Simon/codes/faster_clang

all: executable
executable: amroBuddy3D.c
	$(COMPILER) $(OPTIONS) -o $(EXEC) amroBuddy3D.c -llapack -lblas -lstdc++ -lm
install : all
	cp $(EXEC) $(INSTALLDIR)