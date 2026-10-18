#!/usr/bin/env bash
set -euo pipefail
pipeline_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
pipeline_tools="${pipeline_root}/.tools"
mkdir -p "${pipeline_tools}"
pipeline_expected=60db869a7618783e3ee60183a2b0c6c50ffc67221eabec60800b904991786aaa
pipeline_download="${pipeline_tools}/nextflow.download"
curl --fail --location --retry 3 --output "${pipeline_download}" \
  https://github.com/nextflow-io/nextflow/releases/download/v25.10.4/nextflow-25.10.4-dist
printf '%s  %s\n' "${pipeline_expected}" "${pipeline_download}" | sha256sum --check --status
chmod +x "${pipeline_download}"
mv "${pipeline_download}" "${pipeline_tools}/nextflow"
"${pipeline_tools}/nextflow" -version
printf 'Add this directory to PATH: %s\n' "${pipeline_tools}"
