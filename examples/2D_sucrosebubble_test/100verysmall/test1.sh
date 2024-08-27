#!/usr/bin/env bash





    #>
    #> The MFC prologue prints a summary of the running job and starts a timer.
    #>

    . "/oscar/home/mcarcana/MFC/toolchain/util.sh"

    TABLE_FORMAT_LINE="| * %-14s $MAGENTA%-35s$COLOR_RESET * %-14s $MAGENTA%-35s$COLOR_RESET |\\n"
    TABLE_HEADER="+-----------------------------------------------------------------------------------------------------------+ \\n"
    TABLE_FOOTER="+-----------------------------------------------------------------------------------------------------------+ \\n"
    TABLE_TITLE_FORMAT="| %-105s |\\n"
    TABLE_CONTENT=$(cat <<-END
$(printf "$TABLE_FORMAT_LINE" "Start-time"   "$(date +%T)"                    "Start-date" "$(date +%T)")
$(printf "$TABLE_FORMAT_LINE" "Partition"    "batch"          "Walltime"   "01:00:00")
$(printf "$TABLE_FORMAT_LINE" "Account"      "N/A"          "Nodes"      "1")
$(printf "$TABLE_FORMAT_LINE" "Job Name"     "test1"                        "Engine"     "interactive")
$(printf "$TABLE_FORMAT_LINE" "QoS"          "N/A" "Binary"     "N/A")
$(printf "$TABLE_FORMAT_LINE" "Queue System" "Interactive"                "Email"      "N/A")
END
)

    printf "$TABLE_HEADER"
    printf "$TABLE_TITLE_FORMAT" "MFC case # test1 @ /oscar/home/mcarcana/MFC/examples/2D_sucrosebubble_test/case.py:"
    printf "$TABLE_HEADER"
    printf "$TABLE_CONTENT\\n"
    printf "$TABLE_FOOTER\\n"

        export CUDA_VISIBLE_DEVICES='0'

    t_start=$(date +%s)


ok ":) Loading modules:\n"
cd "/oscar/home/mcarcana/MFC"
. ./mfc.sh load -c o -m c
cd - > /dev/null
echo

    
    ok ":) Running$MAGENTA syscheck$COLOR_RESET:\n"

    if [ 'syscheck' == 'simulation' ]; then
        export CRAY_ACC_MODULE='/oscar/home/mcarcana/MFC/build/staging/7e71bfa6d4/simulation-wg256.lld.exe'
    fi

    cd "/oscar/home/mcarcana/MFC/examples/2D_sucrosebubble_test"

    t_syscheck_start=$(perl -MTime::HiRes=time -e 'printf "%.9f\n", time')


        (set -x;                 mpirun -np 1                                                        "/oscar/home/mcarcana/MFC/build/install/7e71bfa6d4/bin/syscheck")

    
    code=$?

    t_syscheck_stop=$(perl -MTime::HiRes=time -e 'printf "%.9f\n", time')

    if [ $code -ne 0 ]; then
        echo
        error ":( $MAGENTA/oscar/home/mcarcana/MFC/build/install/7e71bfa6d4/bin/syscheck$COLOR_RESET failed with exit code $MAGENTA$code$COLOR_RESET."
        echo
        exit 1
    fi

    unset CRAY_ACC_MODULE



    echo
    
    ok ":) Running$MAGENTA post_process$COLOR_RESET:\n"

    if [ 'post_process' == 'simulation' ]; then
        export CRAY_ACC_MODULE='/oscar/home/mcarcana/MFC/build/staging/88890036a0/simulation-wg256.lld.exe'
    fi

    cd "/oscar/home/mcarcana/MFC/examples/2D_sucrosebubble_test"

    t_post_process_start=$(perl -MTime::HiRes=time -e 'printf "%.9f\n", time')


        (set -x;                 mpirun -np 1                                                        "/oscar/home/mcarcana/MFC/build/install/88890036a0/bin/post_process")

    
    code=$?

    t_post_process_stop=$(perl -MTime::HiRes=time -e 'printf "%.9f\n", time')

    if [ $code -ne 0 ]; then
        echo
        error ":( $MAGENTA/oscar/home/mcarcana/MFC/build/install/88890036a0/bin/post_process$COLOR_RESET failed with exit code $MAGENTA$code$COLOR_RESET."
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
    printf "$TABLE_TITLE_FORMAT" "Finished test1:"
    printf "$TABLE_FORMAT_LINE"  "Total-time:" "$(expr $t_stop - $t_start)s" "Exit Code:" "$code"
    printf "$TABLE_FORMAT_LINE"  "End-time:"   "$(date +%T)"                 "End-date:"  "$(date +%T)"
    printf "$TABLE_FOOTER"

    exit $code

