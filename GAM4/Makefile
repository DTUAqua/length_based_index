.PHONY = all

## Use parallel BLAS for all runs
export OMP_NUM_THREADS=1
export MAX_DURATION=16:00:00
export MAX_MEMORY=8gb

FIRST_LENGTH_GROUP = 1
LAST_LENGTH_GROUP = 33

lengthgroups:= $(shell echo 'cat(formatC($(FIRST_LENGTH_GROUP):$(LAST_LENGTH_GROUP), digits=3, flag="0"))' | R --slave)
resultsfiles:=$(foreach cm,$(lengthgroups),results/out-cm$(cm).RData) 

all: $(resultsfiles) 

COD='Gadus morhua'

## Run the model
RUN_MODEL = R --slave < modelOne.R

MODEL_SRC = modelOne.R

## Targets
results/out-cm%.RData: $(MODEL_SRC)
	mkdir -p results
	export SCRIPT_INPUT="{ SPECIES = $(COD); CMGROUP = $*; OUTFILE='$@' }"; $(RUN_MODEL)


