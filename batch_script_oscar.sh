#!/bin/bash
#SBATCH -N 1
#SBATCH -J run
#SBATCH -o test4.out
#SBATCH -e test4.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --cpu-freq=high
#SBATCH --partition=batch
#SBATCH --mem=2G
#SBATCH --mail-type=all
#SBATCH --mail-user=srijan_neogi@brown.edu

##SBATCH --exclusive
/usr/sbin/ibstat/ | grep Firmware

#Run a command

module load python
module load hpcx-mpi

rm -rf users/sneogi1/MFC/examples/2D_sucrosebubble_test/silo_hdf5/
./mfc.sh run /users/sneogi1/MFC/examples/2D_sucrosebubble_test/case.py -p batch -N 1 -n 1 -g 1 -w 00:20:00 -# test1 -t pre_process -c oscar
./mfc.sh run /users/sneogi1/MFC/examples/2D_sucrosebubble_test/case.py -p batch -N 1 -n 1 -g 1 -w 01:20:00 -# test1 -t simulation -c oscar
./mfc.sh run /users/sneogi1/MFC/examples/2D_sucrosebubble_test/case.py -p batch -N 1 -n 1 -g 1 -w 00:20:00 -# test1 -t post_process -c oscar
