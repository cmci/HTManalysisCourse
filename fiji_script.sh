#!/bin/bash -e
# Read-only filesystem error occurs sometimes
#     Appears to be due to mounting issues or file system limitations
#     Retry if it occurs
BASE_DIR="/g/data/bio-it_centres_course/data/VSVG"
SUB_DIR="${BASE_DIR}/fiji-sub"
LOG_DIR="${BASE_DIR}/fiji-log"
PRE_DIR="${BASE_DIR}/prescreen"
# The first version of the python script gave an error for some plates
#JYTHON_SCRIPT="/g/almf/software/scripts2/measTransportBatch.py"
#JYTHON_SCRIPT="/g/almf/software/scripts2/measTransportBatch2.py"
JYTHON_SCRIPT="/g/almf/software/scripts2/measTransportBatch3.py"

cd ${BASE_DIR}
[ -d fiji-sub ] || mkdir fiji-sub
[ -d fiji-log ] || mkdir fiji-log
[ -d fiji-out ] || mkdir fiji-out

for a in `ls ${BASE_DIR} | grep "\-\-"`
do
    PLATE=${a}

    # Check if prescreen exists
    if [ ! -f ${PRE_DIR}/${PLATE}.csv ]; then
        echo "Skipping plate ${PLATE}: no pre-screen data."
        continue
    fi
    OO="${LOG_DIR}/${PLATE}-out.txt"
    EO="${LOG_DIR}/${PLATE}-err.txt"
    TARGET_DIR="${PLATE}"
    OUTPUT_DIR="fiji-out"

    # Continue if output already exists - no need to re-run
    if [ -f ${OUTPUT_DIR}/${PLATE}--PMall.csv ]; then
        echo "Skipping plate ${PLATE}: results already exist."
        continue
    fi

    # Write the bsub script
    echo "#!/bin/bash
#BSUB -oo \"${OO}\"
#BSUB -eo \"${EO}\"
#BSUB -M 8000
#BSUB -R select[mem>8000] -R rusage[mem=8000]
ulimit -c 0
export PATH=/g/almf/software/bin2:\${PATH}
echo \"The job started:\"
java -version
JYTHON_SCRIPT=${JYTHON_SCRIPT}
TARGET_DIR=${TARGET_DIR}
OUTPUT_DIR=${OUTPUT_DIR}
echo \"Fiji started:\"
fiji --headless --mem=1000m ${JYTHON_SCRIPT} ${TARGET_DIR} ${OUTPUT_DIR}" > ${SUB_DIR}/${PLATE}.bsub

    # Submit the job to the cluster
    chmod +x ${SUB_DIR}/${PLATE}.bsub
    bsub < ${SUB_DIR}/${PLATE}.bsub
done
