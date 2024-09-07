#!/bin/bash

#./mfc.sh run ./examples/3D_hyper_prestress_bubinwater/case.py -p batch -N 1 -n 4 -g 0 -w 01:00:00 -# test1 -t pre_process -c oscar
#./mfc.sh run ./examples/3D_hyper_prestress_bubinwater/case.py -p batch -N 1 -n 4 -g 0 -w 01:00:00 -# test1 -t simulation -c oscar
#./mfc.sh run ./examples/3D_hyper_prestress_bubinwater/case.py -p batch -N 1 -n 4 -g 0 -w 01:00:00 -# test1 -t post_process -c oscar

#./mfc.sh run ./examples/3D_hyper_prestress_bubingel/case.py -p batch -N 1 -n 4 -g 0 -w 01:00:00 -# test1 -t pre_process -c oscar
#./mfc.sh run ./examples/3D_hyper_prestress_bubingel/case.py -p batch -N 1 -n 4 -g 0 -w 01:00:00 -# test1 -t simulation -c oscar
#./mfc.sh run ./examples/3D_hyper_prestress_bubingel/case.py -p batch -N 1 -n 4 -g 0 -w 01:00:00 -# test1 -t post_process -c oscar

./mfc.sh run ../ctr2024/jose-pc/testing/prestress/3d_gelff_prestress_3f.py -p batch -N 1 -n 4 -g 0 -w 01:00:00 -# test1 -t pre_process -c oscar
./mfc.sh run ../ctr2024/jose-pc/testing/prestress/3d_gelff_prestress_3f.py -p batch -N 1 -n 4 -g 0 -w 01:00:00 -# test1 -t simulation -c oscar
./mfc.sh run ../ctr2024/jose-pc/testing/prestress/3d_gelff_prestress_3f.py -p batch -N 1 -n 4 -g 0 -w 01:00:00 -# test1 -t post_process -c oscar

