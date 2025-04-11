#!/bin/bash

patch -p1 -d /opt/bdsim/jairhul-bdsim-c725b23739b8/ < t2k_bdsim_modifications.diff

cd /opt/bdsim/bdsim-build/
cmake ../jairhul-bdsim-c725b23739b8/
make -j8
make install
