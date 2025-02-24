#!/bin/bash

patch -p1 -d /home/bdsim/jairhul-bdsim-c725b23739b8/ < t2k_bdsim_modifications.diff

cd /home/bdsim/bdsim-build/
cmake ../jairhul-bdsim-c725b23739b8/
make -j8
make install
