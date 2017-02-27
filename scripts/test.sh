#!/usr/bin/env bash

set -x  # print commands

./scripts/scanpy.py exdata 
./scripts/scanpy.py examples 
# moignard15, start with rereading
./scripts/scanpy.py moignard15 pca -s 5 --savefigs -r read -v 0
./scripts/scanpy.py moignard15 tsne -s 5 --savefigs -r -v 0
./scripts/scanpy.py moignard15 diffmap -s 5 --savefigs -r  -v 0
./scripts/scanpy.py moignard15 dpt -s 5 --savefigs -r -v 0
./scripts/scanpy.py moignard15 dpt -p basis tsne -s 5 --savefigs -r -v 0
# burczynski06
./scripts/scanpy.py burczynski06 difftest --savefigs -r read -v 0
# paul15pca
./scripts/scanpy.py paul15pca dpt -s 5 --savefigs -r read -v 0
./scripts/scanpy.py paul15pca difftest -o smp dpt_groups names '2,3' -s 5 --savefigs -r -v 0
# toggleswitch
./scripts/scanpy.py toggleswitch sim --savefigs -r read -v 0
./scripts/scanpy.py toggleswitch dpt --savefigs -r -v 0