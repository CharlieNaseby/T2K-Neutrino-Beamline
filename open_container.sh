#!/bin/bash
export DISPLAY=:0
xhost +local:$(whoami)
apptainer shell --bind /tmp/.X11-unix --writable --env DISPLAY=$DISPLAY --env XAUTHORITY=$XAUTHORITY beamline_optimiser/

