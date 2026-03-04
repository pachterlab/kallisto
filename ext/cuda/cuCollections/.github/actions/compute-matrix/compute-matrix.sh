#!/bin/bash
# SPDX-FileCopyrightText: Copyright (c) 2023 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: Apache-2.0
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

set -euo pipefail

write_output() {
  local key="$1"
  local value="$2"
  echo "$key=$value" | tee --append "${GITHUB_OUTPUT:-/dev/null}"
}

explode_std_versions() {
  jq -cr 'map(. as $o | {std: $o.std[]} + del($o.std))'
}

extract_matrix() {
  local file="$1"
  local type="$2"
  local matrix=$(yq -o=json "$file" | jq -cr ".$type")
  write_output "DEVCONTAINER_VERSION" "$(yq -o json "$file" | jq -cr '.devcontainer_version')"
  local nvcc_full_matrix="$(echo "$matrix" | jq -cr '.nvcc' | explode_std_versions )"
  write_output "NVCC_FULL_MATRIX" "$nvcc_full_matrix"
  write_output "CUDA_VERSIONS" "$(echo "$nvcc_full_matrix" | jq -cr '[.[] | .cuda] | unique')"
  write_output "HOST_COMPILERS" "$(echo "$nvcc_full_matrix" | jq -cr '[.[] | .compiler.name] | unique')"
  write_output "PER_CUDA_COMPILER_MATRIX" "$(echo "$nvcc_full_matrix" | jq -cr ' group_by(.cuda + .compiler.name) | map({(.[0].cuda + "-" + .[0].compiler.name): .}) | add')"
}

main() {
  if [ "$1" == "-v" ]; then
    set -x
    shift
  fi

  if [ $# -ne 2 ] || [ "$2" != "pull_request" ]; then
    echo "Usage: $0 [-v] MATRIX_FILE MATRIX_TYPE"
    echo "  -v            : Enable verbose output"
    echo "  MATRIX_FILE   : The path to the matrix file."
    echo "  MATRIX_TYPE   : The desired matrix. Supported values: 'pull_request'"
    exit 1
  fi

  echo "Input matrix file:" >&2
  cat "$1" >&2
  echo "Matrix Type: $2" >&2

  extract_matrix "$1" "$2"
}

main "$@"