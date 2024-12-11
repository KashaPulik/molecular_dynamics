serial: md.cpp
	mpicxx -Wall -o md md.cpp

parallel: pmd.cpp
	mpicxx -Wall -o pmd pmd.cpp

all: serial parallel