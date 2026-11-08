#!/usr/bin/env bash
# Compile the VFIFE core + C shim to a single-file ES module for the web app.
# Requires the Emscripten SDK (emcc); falls back to ~/emsdk if not on PATH.
set -euo pipefail
cd "$(dirname "$0")/.."

if ! command -v em++ >/dev/null 2>&1; then
    if [ -f "$HOME/emsdk/emsdk_env.sh" ]; then
        # shellcheck disable=SC1091
        source "$HOME/emsdk/emsdk_env.sh" >/dev/null 2>&1
    else
        echo "em++ not found; install emsdk or add it to PATH" >&2
        exit 1
    fi
fi

mkdir -p web/src/solver

em++ -O3 -std=c++20 -I include -I headers \
    src/core/particle.cpp src/numerics/quadrature.cpp src/material/material.cpp \
    src/section/section.cpp src/section/grid_section.cpp src/element/element.cpp \
    src/solver/struct_system.cpp \
    wasm/vfife_shim.cpp \
    -s MODULARIZE=1 \
    -s EXPORT_ES6=1 \
    -s ENVIRONMENT=web,worker \
    -s SINGLE_FILE=1 \
    -s ALLOW_MEMORY_GROWTH=1 \
    -s EXPORTED_RUNTIME_METHODS=HEAPF64,HEAP32,UTF8ToString \
    -o web/src/solver/vfife.js

echo "wrote web/src/solver/vfife.js"
