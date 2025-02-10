#!/bin/bash

command="julia --project=. start_MSI_GUI.jl"

case "$OSTYPE" in
  linux*)
    $command
    ;;
  darwin*)
    $command
    ;;
  cygwin*|msys*|win32)
    julia --project=. start_MSI_GUI.jl
    ;;
  *)
    echo "Unsupported OS"
    ;;
esac
