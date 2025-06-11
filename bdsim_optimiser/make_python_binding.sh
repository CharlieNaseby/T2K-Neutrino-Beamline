#!/bin/bash

rm Interface.cpython-39-x86_64-linux-gnu.so
rm -r build/*
python binder.py build_ext --inplace
