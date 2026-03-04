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

# Ensure two arguments are provided
if [ $# -ne 2 ]; then
  echo "Usage: $0 <before-file> <after-file>" >&2
  exit 1
fi

# Print the contents of the before file
echo "=== Contents of $1 ===" >&2
cat $1 >&2
echo "=== End of $1 ===" >&2

# Print the contents of the after file
echo "=== Contents of $2 ==="  >&2
cat $2 >&2
echo "=== End of $2 ===" >&2

# Extract compile requests and cache hits from the before and after files
requests_before=$(awk '/^[ \t]*Compile requests[ \t]+[0-9]+/ {print $3}' "$1")
hits_before=$(awk '/^[ \t]*Cache hits[ \t]+[0-9]+/ {print $3}' "$1")
requests_after=$(awk '/^[ \t]*Compile requests[ \t]+[0-9]+/ {print $3}' "$2")
hits_after=$(awk '/^[ \t]*Cache hits[ \t]+[0-9]+/ {print $3}' "$2")

# Calculate the differences to find out how many new requests and hits
requests_diff=$((requests_after - requests_before))
hits_diff=$((hits_after - hits_before))

echo "New Compile Requests: $requests_diff" >&2
echo "New Hits: $hits_diff" >&2

# Calculate and print the hit rate
if [ $requests_diff -eq 0 ]; then
    echo "No new compile requests, hit rate is not applicable"
else
    hit_rate=$(awk -v hits=$hits_diff -v requests=$requests_diff 'BEGIN {printf "%.2f", hits/requests * 100}')
    echo "sccache hit rate: $hit_rate%" >&2
    echo "$hit_rate"
fi