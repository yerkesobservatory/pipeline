#/bin/bash

# PipeDailyRun
# ============
#
# Script to automatically run pipeline services. The following things are run:
#   * Make master bias / darks / flats
#   * Run all of today's data by using PipeExecuteAutoDay

### Setup
export PYTHONPATH=/data/scripts/DataReduction/source

### Run Masters
cd /data/scripts/DataReduction
/usr/local/bin/python source/drp/pipeline.py --loglevel DEBUG --logfile PipeLineLog.txt --pipemode masterbias pipeconf_stonedge_auto.txt >> AstroLog.txt 2>&1
/usr/local/bin/python source/drp/pipeline.py --loglevel DEBUG --logfile PipeLineLog.txt --pipemode masterdark pipeconf_stonedge_auto.txt >> AstroLog.txt 2>&1
/usr/local/bin/python source/drp/pipeline.py --loglevel DEBUG --logfile PipeLineLog.txt --pipemode masterflat pipeconf_stonedge_auto.txt >> AstroLog.txt 2>&1

### Run Pipeline
./PipeExecuteAutoDay.py >> AstroLog.txt 2>&1