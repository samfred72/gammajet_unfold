#!/bin/bash
# Build the multijet balance analysis (standalone executable, not part of libgammajet_unfold).
cd "$(dirname "$0")"
g++ analysis.cc style.cc -std=c++17 `root-config --cflags --libs` -o analysis
