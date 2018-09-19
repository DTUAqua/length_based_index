.PHONY = all

## Use parallel BLAS for all runs
export OMP_NUM_THREADS=2
export MAX_DURATION=01:00:00:00
export MAX_MEMORY=35gb
export QUEUE=hpc

FIRST_LENGTH_GROUP = 1
LAST_LENGTH_GROUP = 38

lengthgroups:= $(shell echo 'cat(formatC($(FIRST_LENGTH_GROUP):$(LAST_LENGTH_GROUP), digits=3, flag="0"))' | R --slave)
resultsQ14files:=$(foreach cm,$(lengthgroups),results/cod-Q14-cm$(cm).RData)
animQ1files:=$(foreach cm,$(lengthgroups),anim/cod-Q1-cm$(cm).pdf)
animQ4files:=$(foreach cm,$(lengthgroups),anim/cod-Q4-cm$(cm).pdf)
postQ14files:=$(foreach cm,$(lengthgroups),post/cod-Q14-cm$(cm).RData)

all: $(resultsQ14files)

all_anim: $(animQ1files) $(animQ4files)

all_post: $(postQ14files)

join:
	R --slave < join.R

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

anim/cod-Q1-cm%.pdf: $(MODEL_BIN)
	mkdir -p anim
	export SCRIPT_INPUT="{ QUARTERS='1'; RESFILE='results/cod-Q14-cm$*.RData'; OUTFILE='$@' }"; R --slave < animation.R

anim/cod-Q4-cm%.pdf: $(MODEL_BIN)
	mkdir -p anim
	export SCRIPT_INPUT="{ QUARTERS='4'; RESFILE='results/cod-Q14-cm$*.RData'; OUTFILE='$@' }"; R --slave < animation.R

post/cod-Q14-cm%.RData: $(MODEL_BIN)
	mkdir -p post
	export SCRIPT_INPUT="{ RESFILE='results/cod-Q14-cm$*.RData'; OUTFILE='$@' }"; R --slave < postprocess.R

## Optional: Send mail when done + link to results
publish: results
	echo "Copy result files here..." > ~/public_html/results.txt
	chmod 755 ~/public_html/results.txt
	echo "Results are at www.student.dtu.dk/~"$$USER | mail -s "Results have been updated" $$USER@dtu.dk

wget-results:
	wget http://www.student.dtu.dk/~kaskr/EBcod_results.zip
	unzip EBcod_results.zip
