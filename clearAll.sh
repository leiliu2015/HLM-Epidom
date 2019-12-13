#!/bin/bash

fds=(A-23 I-14 R-11 Fig-S6 Fig-S7 Fig-S8 Fig-S9)
aps=(rgs rho srf asp irg)

for ((idx=0; "$idx" < 7; idx++))
do
    fd=${fds[$idx]}
    echo ${fd}
    #
    cd ./${fd}
        cd ./MD
            if [ -f gnm.t0.sys ]
            then
                rm ./gnm.t0.*
            fi
        cd ../
        #
        if [ -f gnm.t0-0.rc*.cm ]
        then
            rm ./gnm.t0-0.rc*.cm
        fi
        #
        for ((kdx=0; "$kdx" < 5; kdx++))
        do
            apx=${aps[$kdx]}
            if [ -f gnm.t0-0.2016zhuang418.${apx} ]
            then
                rm ./gnm.t0-0.2016zhuang418.${apx}
            fi
        done
    cd ../
done

