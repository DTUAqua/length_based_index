.PHONY = all

## Use parallel BLAS for all runs
export OMP_NUM_THREADS=8
export MAX_DURATION=01:00:00:00
export MAX_MEMORY=35gb
export QUEUE=computebigmem

FIRST_LENGTH_GROUP = 5
LAST_LENGTH_GROUP = 100

lengthgroups:= $(shell echo 'cat(formatC($(FIRST_LENGTH_GROUP):$(LAST_LENGTH_GROUP), digits=3, flag="0"))' | R --slave)
resultsQ14files:=$(foreach cm,$(lengthgroups),results/cod-Q14-cm$(cm).RData)

all: $(resultsQ14files)

COD='Gadus morhua'

## Run the model
RUN_MODEL = R --slave < model.R

## Build the model
MODEL_SRC = model.cpp
MODEL_BIN = model.so
$(MODEL_BIN) : $(MODEL_SRC)
	echo "TMB:::compile('$(MODEL_SRC)', '-Ofast')" | R --slave

## Q14 Targets
results/cod-Q14-cm%.RData: $(MODEL_BIN)
	mkdir -p results
	export SCRIPT_INPUT="{ SPECIES = $(COD); CMGROUP = $*; OUTFILE='$@' }"; $(RUN_MODEL)

## Optional: Send mail when done + link to results
publish: results
	echo "Copy result files here..." > ~/public_html/results.txt
	chmod 755 ~/public_html/results.txt
	echo "Results are at www.student.dtu.dk/~"$$USER | mail -s "Results have been updated" $$USER@dtu.dk
