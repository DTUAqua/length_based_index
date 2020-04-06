.PHONY = all

## Use parallel BLAS for all runs
export OMP_NUM_THREADS=1
export MAX_DURATION=24:00:00
export MAX_MEMORY=6gb

FIRST_LENGTH_GROUP = 1
LAST_LENGTH_GROUP = 9

FIRST_SETTING = 1
LAST_SETTING = 20

lengthgroups:= $(shell echo 'cat(formatC($(FIRST_LENGTH_GROUP):$(LAST_LENGTH_GROUP), digits=3, flag="0"))' | R --slave)
resultsfiles:=$(foreach cm,$(lengthgroups),alldenssplit-$(cm).RData)
settings:= $(shell echo 'cat(formatC($(FIRST_SETTING):$(LAST_SETTING), digits=3, flag="0"))' | R --slave)
settingfiles:=$(foreach cm,$(settings),ageindices_duchon-$(cm).RData)

all: $(resultsfiles) 

ageindices: $(settingfiles)

## Run the model
RUN_MODEL = R --slave < splitalldens.R

MODEL_SRC = splitalldens.R

SETTING_SRC = script_duchonALK-one.R

RUN_SETTINGS = R --slave < script_duchonALK-one.R

## Targets
alldenssplit-%.RData: $(MODEL_SRC)
	export SCRIPT_INPUT="{ MODELNR = $*; OUTFILE='$@' }"; $(RUN_MODEL)

get-results:
	scp cabe@transfer.gbar.dtu.dk:length_based_index/WBcod/*.RData .

ageindices_duchon-%.RData: $(SETTING_SRC)
	export SCRIPT_INPUT="{ SETTING = $*; OUTFILE='$@' }"; $(RUN_SETTINGS)
