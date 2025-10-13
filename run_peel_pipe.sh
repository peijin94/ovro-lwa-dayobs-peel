#!/bin/bash
# Run LWA peeling pipeline in container for wideband data

# Define container paths (these will be mounted)
SOFT_DIR="/fast/peijinz/day-peel/ovro-lwa-dayobs-peel"
DATA_DIR="/fast/peijinz/day-peel/data1"

FRQ_LIST=(59 64 69 73)
DATETIME="20251009_212514"
BCAL_PREFIX="20250920_041508"

# Run the peeling pipeline in container
podman run --rm -it \
  -v ${SOFT_DIR}:/lwasoft:ro \
  -v ${DATA_DIR}:/data:rw \
  -w /data \
  peijin/lwa-solar-pipehost:v202510 \
  python3 /lwasoft/peel_pipe.py \
    --ms-list \
      /data/slow/${DATETIME}_${FRQ_LIST[0]}MHz.ms \
      /data/slow/${DATETIME}_${FRQ_LIST[1]}MHz.ms \
      /data/slow/${DATETIME}_${FRQ_LIST[2]}MHz.ms \
      /data/slow/${DATETIME}_${FRQ_LIST[3]}MHz.ms \
    --bcal-list \
      /data/bcal/${BCAL_PREFIX}_${FRQ_LIST[0]}MHz.bcal \
      /data/bcal/${BCAL_PREFIX}_${FRQ_LIST[1]}MHz.bcal \
      /data/bcal/${BCAL_PREFIX}_${FRQ_LIST[2]}MHz.bcal \
      /data/bcal/${BCAL_PREFIX}_${FRQ_LIST[3]}MHz.bcal \
    --output-prefix "peel_${DATETIME}" \
    --output-dir "/data/slow" \
  > peel_proc.log

echo "Peeling pipeline completed! Check peel_proc.log for details."

