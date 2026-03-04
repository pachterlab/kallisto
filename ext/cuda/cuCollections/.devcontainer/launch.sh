#! /usr/bin/env bash
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

launch_devcontainer() {

    # Ensure we're in the repo root
    cd "$( cd "$( dirname "$(realpath -m "${BASH_SOURCE[0]}")" )" && pwd )/..";

    if [[ -z $1 ]] || [[ -z $2 ]]; then
        echo "Usage: $0 [CUDA version] [Host compiler]"
        echo "Example: $0 12.1 gcc12"
        return 1
    fi

    local cuda_version="$1"
    local host_compiler="$2"
    local workspace="$(basename "$(pwd)")";
    local tmpdir="$(mktemp -d)/${workspace}";
    local path="$(pwd)/.devcontainer/cuda${cuda_version}-${host_compiler}";

    mkdir -p "${tmpdir}";
    mkdir -p "${tmpdir}/.devcontainer";
    cp -arL "$path/devcontainer.json" "${tmpdir}/.devcontainer";
    sed -i "s@\${localWorkspaceFolder}@$(pwd)@g" "${tmpdir}/.devcontainer/devcontainer.json";
    path="${tmpdir}";

    local hash="$(echo -n "${path}" | xxd -pu - | tr -d '[:space:]')";
    local url="vscode://vscode-remote/dev-container+${hash}/home/coder/cuCollections";

    echo "devcontainer URL: ${url}";

    local launch="";
    if type open >/dev/null 2>&1; then
        launch="open";
    elif type xdg-open >/dev/null 2>&1; then
        launch="xdg-open";
    fi

    if [ -n "${launch}" ]; then
        code --new-window "${tmpdir}";
        exec "${launch}" "${url}" >/dev/null 2>&1;
    fi
}

launch_devcontainer "$@";