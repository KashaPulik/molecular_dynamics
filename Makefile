SERIAL_FILES=./src/initialization.cpp ./src/md.cpp ./src/physical.cpp ./src/serial_model.cpp
PARALLEL_FILES=./src/initialization.cpp ./src/pmd.cpp ./src/physical.cpp ./src/parallel_model.cpp

all: serial parallel

serial: $(SERIAL_FILES)
	mpicxx -Wall -o md $(SERIAL_FILES)

parallel: $(PARALLEL_FILES)
	mpicxx -Wall -o pmd $(PARALLEL_FILES)

clean:
	rm -rf md pmd slurm*