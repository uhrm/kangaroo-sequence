#!/bin/sh
# Install the tools used by this project into ~/.local/bin and set up the Python environment.
set -eu

UV_VERSION=0.12.19
TECTONIC_VERSION=0.17.0

mkdir -p ~/.local/bin

# uv manages Python and the project's virtual environment.
curl -LsSf "https://astral.sh/uv/$UV_VERSION/install.sh" | sh
~/.local/bin/uv sync

# Tectonic builds the report. The statically linked (musl) build avoids missing system libraries.
curl -fsSL "https://github.com/tectonic-typesetting/tectonic/releases/download/tectonic%40$TECTONIC_VERSION/tectonic-$TECTONIC_VERSION-$(uname -m)-unknown-linux-musl.tar.gz" \
    | tar xz -C ~/.local/bin
~/.local/bin/tectonic --version
