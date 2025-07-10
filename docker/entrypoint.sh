#!/bin/bash
set -e

if [[ "$1" == "rstudio" ]]; then
    echo "Starting RStudio Server..."
    exec /init
    # exec /usr/lib/rstudio-server/bin/rserver --server-daemonize=0
else
    echo "Running arbitrary command: $@"
    exec "$@"
fi