OPTIONS = -O2 -Wall
EXEC = amroBuddy
#INSTALLDIR := $(HOME)/bin/
#COMPILER = gcc
COMPILER = /Users/Simon/codes/faster_clang

all: executable
executable: amroBuddy.c
	$(COMPILER) $(OPTIONS) -o $(EXEC) amroBuddy.c -llapack -lblas -lstdc++ -lm
fast: amroFastBuddy.c
	$(COMPILER) $(OPTIONS) -o fastBuddy amroFastBuddy.c -llapack -lblas -lstdc++ -lm

# install : all
# 	cp $(EXEC) $(INSTALLDIR)