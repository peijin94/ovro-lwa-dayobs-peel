#!/bin/bash
# Run LWA solar pipeline in container

podman run --rm -it \
  -v /fast/peijinz/day-peel/ovro-lwa-dayobs-peel:/lwasoft:ro \
  -v /fast/peijinz/day-peel/data:/data:rw \
  -w /data \
  peijin/lwa-solar-pipehost:v202510 \
  python3 /lwasoft/peel_pipe.py \
    /data/20251010_195103_69MHz.ms \
    /data/20250920_041508_69MHz.bcal \
  > proc.log