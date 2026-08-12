#/bin/bash

# PipeDailyRun
# ============
#
# Script to automatically run pipeline services. The following things are run:
#   * Make master bias / darks / flats
#   * Run all of today's data by using PipeExecuteAutoDay

### Setup
export PATH=/usr/lib64/qt-3.3/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:$HOME/bin:$PATH
# Have extra /darepype path to use non-published version of darepype
export PYTHONPATH=/data/scripts/pipeline/source

PythonVenv=/data/scripts/pipeline/.venv
PipelineDir= /data/scripts/pipeline

DRPath=$PythonVenv/lib/python3.12/site-packages/darepype

### Clear Temp Files
find /data/images/fitsview/temp/images -type f -atime +2 -delete
find /tmp -type f \( ! -user 0 \) -atime +2 -delete

### Run Sort Obs
cd "$PipelineDir"

source "$PythonVenv/bin/activate"

date >> AstroLog.txt

python3 $DRPath/drp/pipeline.py --loglevel DEBUG --logfile PipeLineLog.txt --pipemode sortobs -c config/dconf_stars.txt config/pipeconf_SEO.txt >> AstroLog.txt 2>&1
### Run Masters
python3 $DRPath/drp/pipeline.py --loglevel DEBUG --logfile PipeLineLog.txt --pipemode masterbias -c config/dconf_stars.txt config/pipeconf_SEO.txt >> AstroLog.txt 2>&1
python3 $DRPath/drp/pipeline.py --loglevel DEBUG --logfile PipeLineLog.txt --pipemode masterdark -c config/dconf_stars.txt config/pipeconf_SEO.txt >> AstroLog.txt 2>&1
python3 $DRPath/drp/pipeline.py --loglevel DEBUG --logfile PipeLineLog.txt --pipemode masterflat -c config/dconf_stars.txt config/pipeconf_SEO.txt >> AstroLog.txt 2>&1

### Run Pipeline
python3 ./PipeExecuteAutoDay.py >> AstroLog.txt 2>&1

### Copy queue data
python3 /usr/local/lib/python3.6/site-packages/darepype/drp/pipeline.py --loglevel DEBUG --logfile PipeLineLog.txt --pipemode copyqueue -c config/dconf_stars.txt config/pipeconf_SEO.txt >> AstroLog.txt 2>&1
### Run the queue pipeline
for dpr in $(ls /data/images/queue/A_Test/piperuns/*`date "+_%Y-%m-%d_"`*.txt) ; do darepyperun.py $dpr ; done >> /data/scripts/pipeline/QueueLog.txt 2>&1