#!/bin/bash
#./mfc.sh run /scratch/bciv/mcarcanabarbosa/ctr2024/finalruns/dfd2024/liq/nopc/case.py -p gpuA40x4 -N 1 -n 2 -g 2 -w 02:00:00 -# liq25 -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/mcarcanabarbosa/ctr2024/finalruns/dfd2024/liq/nopc/case.py -p gpuA40x4 -N 1 -n 2 -g 2 -w 18:00:00 -# liq25 -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/mcarcanabarbosa/ctr2024/finalruns/dfd2024/liq/nopc/case.py -p gpuA40x4 -N 1 -n 2 -g 2 -w 01:00:00 -# liq25 -t post_process -a bciv-delta-gpu -c delta

#./mfc.sh run /scratch/bciv/mcarcanabarbosa/ctr2024/finalruns/dfd2024/gel/case.py -e batch -p gpuA40x4 -N 1 -n 4 -g 4 -w 01:00:00 -# liq25 -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/mcarcanabarbosa/ctr2024/finalruns/dfd2024/gel/case.py -e batch -p gpuA40x4 -N 1 -n 4 -g 4 -w 01:00:00 -# liq25 -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/mcarcanabarbosa/ctr2024/finalruns/dfd2024/gel/case.py -p gpuA40x4 -N 1 -n 4 -g 4 -w 01:00:00 -# liq25 -t post_process -a bciv-delta-gpu -c delta

./mfc.sh run /scratch/bciv/mcarcanabarbosa/dfd2024/pc_hypo_wall/3d_pc_3f_lff_half.py -p gpuA40x4 -N 1 -n 2 -g 2 -w 02:00:00 -# pclff -t pre_process -a bciv-delta-gpu -c delta
./mfc.sh run /scratch/bciv/mcarcanabarbosa/dfd2024/pc_hypo_wall/3d_pc_3f_lff_half.py -p gpuA40x4 -N 1 -n 2 -g 2 -w 02:00:00 -# pclff -t simulation -a bciv-delta-gpu -c delta
./mfc.sh run /scratch/bciv/mcarcanabarbosa/dfd2024/pc_hypo_wall/3d_pc_3f_lff_half.py -p gpuA40x4 -N 1 -n 2 -g 2 -w 02:00:00 -# pclff -t post_process -a bciv-delta-gpu -c delta

#./mfc.sh run ./examples/1D_hyper_impact/case.py -p gpuA40x4 -N 1 -n 1 -g 1 -w 01:00:00 -# erweak1 -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run ./examples/1D_hyper_impact/case.py -p gpuA40x4 -N 1 -n 1 -g 1 -w 01:00:00 -# erweak2 -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run ./examples/1D_hyper_impact/case.py -p gpuA40x4 -N 1 -n 1 -g 1 -w 01:00:00 -# erweak3 -t post_process -a bciv-delta-gpu -c delta
#./mfc.sh run ./hypertest/1D_hypo_impact_weak/case.py -p gpuA40x4 -N 1 -n 1 -g 1 -w 01:00:00 -# oweak1 -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run ./hypertest/1D_hypo_impact_weak/case.py -p gpuA40x4 -N 1 -n 1 -g 1 -w 01:00:00 -# oweak2 -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run ./hypertest/1D_hypo_impact_weak/case.py -p gpuA40x4 -N 1 -n 1 -g 1 -w 01:00:00 -# oweak3 -t post_process -a bciv-delta-gpu -c delta


#./mfc.sh run ./examples/3D_hyper_prestress_bubinwater/case.py -p batch -N 1 -n 1 -g 0 -w 01:00:00 -# test1 -t pre_process -c oscar
#./mfc.sh run ./examples/3D_hyper_prestress_bubinwater/case.py -p batch -N 1 -n 1 -g 0 -w 01:00:00 -# test1 -t simulation -c oscar
#./mfc.sh run ./examples/3D_hyper_prestress_bubinwater/case.py -p batch -N 1 -n 1 -g 0 -w 01:00:00 -# test1 -t post_process -c oscar

##./mfc.sh run ./examples/3D_hyper_prestress_bubingel/case.py -p batch -N 1 -n 1 -g 0 -w 01:00:00 -# test1 -t pre_process -c oscar
#./mfc.sh run ./examples/3D_hyper_prestress_bubingel/case.py -p batch -N 1 -n 1 -g 0 -w 01:00:00 -# test1 -t simulation -c oscar
##./mfc.sh run ./examples/3D_hyper_prestress_bubingel/case.py -p batch -N 1 -n 1 -g 0 -w 01:00:00 -# test1 -t post_process -c oscar

#./mfc.sh run ../ctr2024/jose-pc/testing/prestress/3d_gelff_prestress_3f.py -p batch -N 1 -n 4 -g 0 -w 01:00:00 -# test1 -t pre_process -c oscar
#./mfc.sh run ../ctr2024/jose-pc/testing/prestress/3d_gelff_prestress_3f.py -p batch -N 1 -n 4 -g 0 -w 01:00:00 -# test1 -t simulation -c oscar
#./mfc.sh run ../ctr2024/jose-pc/testing/prestress/3d_gelff_prestress_3f.py -p batch -N 1 -n 4 -g 0 -w 01:00:00 -# test1 -t post_process -c oscar

#./mfc.sh run /scratch/bciv/mcarcanabarbosa/ctr2024/finalruns/elcom/afterctr/prestress/3d_gelff_prestress_3f.py -e batch -p gpuA40x4 -N 1 -n 1 -g 0 -w 01:00:00 -# test1 -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/mcarcanabarbosa/ctr2024/finalruns/elcom/afterctr/prestress/3d_gelff_prestress_3f.py -e batch -p gpuA40x4 -N 1 -n 1 -g 0 -w 01:00:00 -# test1 -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/mcarcanabarbosa/ctr2024/finalruns/elcom/afterctr/prestress/3d_gelff_prestress_3f.py -e batch -p gpuA40x4 -N 1 -n 1 -g 0 -w 01:00:00 -# test1 -t post_process -a bciv-delta-gpu -c delta
