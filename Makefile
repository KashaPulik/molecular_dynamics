serial: md.cpp
	mpicxx -Wall -o md md.cpp

parallel: pmd.cpp
	mpicxx -Wall -o pmd pmd.cpp

research: mdm.cpp
	mpicxx -Wall -o mdm mdm.cpp

all: serial parallel