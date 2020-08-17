#/bin/bash

# PipeDailyRun
# ============
#
# Script to automatically run pipeline services. The following things are run:
#   * Make master bias / darks / flats
#   * Run all of today's data by using PipeExecuteAutoDay

### Setup
export PATH=/usr/lib64/qt-3.3/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:$HOME/bin:$PATH
export PYTHONPATH=/data/scripts/DataReduction/source

DRPath=/usr/local/lib/python3.6/site-packages/darepype

### Run Masters
cd /data/scripts/DataReduction
/usr/local/bin/python3 $DRPath/drp/pipeline.py --loglevel DEBUG --logfile PipeLineLog.txt --pipemode masterbias pipeconf_stonedge_auto.txt >> AstroLog.txt 2>&1
/usr/local/bin/python3 $DRPath/drp/pipeline.py --loglevel DEBUG --logfile PipeLineLog.txt --pipemode masterdark pipeconf_stonedge_auto.txt >> AstroLog.txt 2>&1
/usr/local/bin/python3 $DRPath/drp/pipeline.py --loglevel DEBUG --logfile PipeLineLog.txt --pipemode masterflat pipeconf_stonedge_auto.txt >> AstroLog.txt 2>&1

### Run Pipeline
./PipeExecuteAutoDay.py >> AstroLog.txt 2>&1
