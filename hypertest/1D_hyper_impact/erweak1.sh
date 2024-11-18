#!/usr/bin/env bash





    #>
    #> The MFC prologue prints a summary of the running job and starts a timer.
    #>

    . "/u/mcarcanabarbosa/MFC/toolchain/util.sh"

    TABLE_FORMAT_LINE="| * %-14s $MAGENTA%-35s$COLOR_RESET * %-14s $MAGENTA%-35s$COLOR_RESET |\\n"
    TABLE_HEADER="+-----------------------------------------------------------------------------------------------------------+ \\n"
    TABLE_FOOTER="+-----------------------------------------------------------------------------------------------------------+ \\n"
    TABLE_TITLE_FORMAT="| %-105s |\\n"
    TABLE_CONTENT=$(cat <<-END
$(printf "$TABLE_FORMAT_LINE" "Start-time"   "$(date +%T)"                    "Start-date" "$(date +%T)")
$(printf "$TABLE_FORMAT_LINE" "Partition"    "gpuA40x4"          "Walltime"   "01:00:00")
$(printf "$TABLE_FORMAT_LINE" "Account"      "bciv-delta-gpu"          "Nodes"      "1")
$(printf "$TABLE_FORMAT_LINE" "Job Name"     "erweak1"                        "Engine"     "interactive")
$(printf "$TABLE_FORMAT_LINE" "QoS"          "N/A" "Binary"     "N/A")
$(printf "$TABLE_FORMAT_LINE" "Queue System" "Interactive"                "Email"      "N/A")
END
)

    printf "$TABLE_HEADER"
    printf "$TABLE_TITLE_FORMAT" "MFC case # erweak1 @ /u/mcarcanabarbosa/MFC/examples/1D_hyper_impact/case.py:"
    printf "$TABLE_HEADER"
    printf "$TABLE_CONTENT\\n"
    printf "$TABLE_FOOTER\\n"

        export CUDA_VISIBLE_DEVICES='1'

    t_start=$(date +%s)


ok ":) Loading modules:\n"
cd "/u/mcarcanabarbosa/MFC"
. ./mfc.sh load -c d -m g
cd - > /dev/null
echo

# Fixes Delta not being able to find core library file
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/spack/deltas11-2023-03/apps/linux-rhel8-zen3/nvhpc-22.11/openmpi-4.1.5-nzb4n4r/lib/

    
    ok ":) Running$MAGENTA syscheck$COLOR_RESET:\n"

    if [ 'syscheck' == 'simulation' ]; then
        export CRAY_ACC_MODULE='/u/mcarcanabarbosa/MFC/build/staging/ea58b6820d/simulation-wg256.lld.exe'
    fi

    cd "/u/mcarcanabarbosa/MFC/examples/1D_hyper_impact"

    t_syscheck_start=$(perl -MTime::HiRes=time -e 'printf "%.9f\n", time')


        (set -x;                                                mpirun -np 1                                                        "/u/mcarcanabarbosa/MFC/build/install/ea58b6820d/bin/syscheck")

    
    code=$?

    t_syscheck_stop=$(perl -MTime::HiRes=time -e 'printf "%.9f\n", time')

    if [ $code -ne 0 ]; then
        echo
        error ":( $MAGENTA/u/mcarcanabarbosa/MFC/build/install/ea58b6820d/bin/syscheck$COLOR_RESET failed with exit code $MAGENTA$code$COLOR_RESET."
        echo
        exit 1
    fi

    unset CRAY_ACC_MODULE



    echo
    
    ok ":) Running$MAGENTA pre_process$COLOR_RESET:\n"

    if [ 'pre_process' == 'simulation' ]; then
        export CRAY_ACC_MODULE='/u/mcarcanabarbosa/MFC/build/staging/25eb1bf60b/simulation-wg256.lld.exe'
    fi

    cd "/u/mcarcanabarbosa/MFC/examples/1D_hyper_impact"

    t_pre_process_start=$(perl -MTime::HiRes=time -e 'printf "%.9f\n", time')


        (set -x;                                                mpirun -np 1                                                        "/u/mcarcanabarbosa/MFC/build/install/25eb1bf60b/bin/pre_process")

    
    code=$?

    t_pre_process_stop=$(perl -MTime::HiRes=time -e 'printf "%.9f\n", time')

    if [ $code -ne 0 ]; then
        echo
        error ":( $MAGENTA/u/mcarcanabarbosa/MFC/build/install/25eb1bf60b/bin/pre_process$COLOR_RESET failed with exit code $MAGENTA$code$COLOR_RESET."
        echo
        exit 1
    fi

    unset CRAY_ACC_MODULE



    echo


    #>
    #> The MFC epilogue stops the timer and prints the execution summary. It also
    #> performs some cleanup and housekeeping tasks before exiting.
    #>

    code=$?

    t_stop="$(date +%s)"

    printf "$TABLE_HEADER"
    printf "$TABLE_TITLE_FORMAT" "Finished erweak1:"
    printf "$TABLE_FORMAT_LINE"  "Total-time:" "$(expr $t_stop - $t_start)s" "Exit Code:" "$code"
    printf "$TABLE_FORMAT_LINE"  "End-time:"   "$(date +%T)"                 "End-date:"  "$(date +%T)"
    printf "$TABLE_FOOTER"

    exit $code

