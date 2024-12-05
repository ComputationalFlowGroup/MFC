#./mfc.sh run ./examples/2D_sucrose_plasticity/case.py -p gpuA40x4 -N 1 -n 1 -g 1 -w 01:00:00 -# test1 -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run ./examples/2D_sucrose_plasticity/case.py -p gpuA40x4 -N 1 -n 1 -g 1 -w 01:00:00 -# test1 -t simulation -a bciv-delta-gpu -c delta
./mfc.sh run ./examples/2D_sucrose_plasticity/case.py -p gpuA40x4 -N 1 -n 1 -g 1 -w 01:00:00 -# test1 -t post_process -a bciv-delta-gpu -c delta
