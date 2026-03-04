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

# This script prints the sccache hit rate between two calls to sccache --show-stats.
# It should be sourced in your script before and after the operations you want to profile,
# with the 'start' or 'end' argument respectively.

mode=$1

if [[ "$mode" != "start" && "$mode" != "end" ]]; then
    echo "Invalid mode: $mode"
    echo "Usage: $0 {start|end}"
    exit 1
fi

case $mode in
  start)
    export SCCACHE_START_HITS=$(sccache --show-stats | awk '/^[ \t]*Cache hits[ \t]+[0-9]+/ {print $3}')
    export SCCACHE_START_MISSES=$(sccache --show-stats | awk '/^[ \t]*Cache misses[ \t]+[0-9]+/ {print $3}')
    ;;
  end)
    if [[ -z ${SCCACHE_START_HITS+x} || -z ${SCCACHE_START_MISSES+x} ]]; then
        echo "Error: start stats not collected. Did you call this script with 'start' before your operations?"
        exit 1
    fi

    final_hits=$(sccache --show-stats | awk '/^[ \t]*Cache hits[ \t]+[0-9]+/ {print $3}')
    final_misses=$(sccache --show-stats | awk '/^[ \t]*Cache misses[ \t]+[0-9]+/ {print $3}')
    hits=$((final_hits - SCCACHE_START_HITS))
    misses=$((final_misses - SCCACHE_START_MISSES))
    total=$((hits + misses))

    prefix=""
    if [ ${GITHUB_ACTIONS:-false} = "true" ]; then
      prefix="::notice::"
    fi

    if (( total > 0 )); then
      hit_rate=$(awk -v hits="$hits" -v total="$total" 'BEGIN { printf "%.2f", (hits / total) * 100 }')
      echo ${prefix}"sccache hits: $hits | misses: $misses | hit rate: $hit_rate%"
    else
      echo ${prefix}"sccache stats: N/A No new compilation requests"
    fi
    unset SCCACHE_START_HITS
    unset SCCACHE_START_MISSES
    ;;
esac