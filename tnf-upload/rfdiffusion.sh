#!/usr/bin/env bash
# RFdiffusion via the maintained Docker image.
#
# The image's ENTRYPOINT already invokes run_inference.py, so pass ONLY Hydra
# overrides -- no `python`, no script path, or Hydra tries to parse them as
# config keys. Do not override the container's workdir either; the entrypoint
# resolves its script path relative to /app/RFdiffusion.
#
# $PWD is mounted at /work, so pass container-absolute paths: /work/target/...
#
# Two dirs must be mounted writable because we run as the host UID and /app is
# root-owned in the image:
#   /app/RFdiffusion/schedules  - IGSO3 diffusion schedule cache, mkdir'd on
#                                 first run (model_runners.py:129). Without
#                                 this the run dies with PermissionError before
#                                 any sampling. Keep it across runs: it saves
#                                 several minutes of "Calculating IGSO3".
#   hydra.run.dir               - pass =/work/.hydra as an override.
set -euo pipefail

mkdir -p "$PWD/schedules"

exec docker run --rm --gpus all \
  --user "$(id -u):$(id -g)" \
  -v "$PWD:/work" \
  -v "$PWD/schedules:/app/RFdiffusion/schedules" \
  -v "$HOME/tnf-track3/RFdiffusion/models:/app/RFdiffusion/models:ro" \
  rosettacommons/rfdiffusion:latest "$@"
