# read in {k_ij}.
proc fread {kapFile} {
    # 
    set fp [open $kapFile r]
    fconfigure $fp -buffering line
    # 
    set kap [list]
    gets $fp data
    while {![eof $fp]} {
        lappend kap $data
        gets $fp data
    }
    # 
    close $fp
    return $kap
}

# starting time 
set st [clock format [clock seconds] -format %Y,%b.%d-%H:%M:%S] 

# set random seed for tcl (srand) and espresso (t_random, bit_random) 
set rSeed 12345 
expr srand($rSeed) 

set ncore 1
set rss {}
for {set i 0} {$i < $ncore} {incr i} {  
    lappend rss $rSeed
} 
eval [format "t_random   seed %s" $rss] 
eval [format "bit_random seed %s" $rss] 

# sys params 
set num 100
set bod 1.0
set ksi 1.0
set tep 1.0 
set box [expr max(15.0, int(pow($num, 0.6)) + 6.0)]
set cob [expr $box/2.0] 
# sim params
setmd box_l $box $box $box 
setmd periodic 1 1 1
setmd skin 0.4 
setmd time_step 0.01 
# job
set tdx 0
set jobName [format "gnm.t%d" $tdx]
set kapName hlm.km

####### Set up system ##########################################  
# initial cfg 
source hlm.iniCfg

### interactions 
# backbone bonds
set dr 0.0
inter 0 harmonic 3.0 0.0
for {set i 1} {$i < $num} {incr i} {
    part $i bond 0 [expr $i-1]
}

# loop bonds
set km [fread $kapName]
set ksf 1.0 
set kdx 1
#
for {set i 0} {$i < $num} {incr i} {
    set ip2 [expr $i+2]
    for {set j $ip2} {$j < $num} {incr j} {
        set kb [lindex $km $i $j]
        #
        if {$kb > 0} {
            inter $kdx harmonic [expr $kb*$ksf] $dr
            part $i bond $kdx $j
            #
            incr kdx
        }
    }
}

# nonbonded interactions
set esf 0.35
inter 0 0 lj-gen [expr 8.69692e-01*$esf] 1.0 2.5 0.008175222784 0.0 12 6 1.0 2.0
inter 0 1 lj-gen [expr 8.13336e-01*$esf] 1.0 2.5 0.008175222784 0.0 12 6 1.0 2.0
inter 1 1 lj-gen [expr 1.25405e+00*$esf] 1.0 2.5 0.008175222784 0.0 12 6 1.0 2.0

####### Set up simulation ########################################## 
# thermostat 
set langevin_temperature $tep  
set langevin_friction $ksi 
thermostat langevin $langevin_temperature $langevin_friction 
integrate set nvt 

# -----------------------------------------------------------
# production run (totalEqulibriumSteps = equ_rounds * equ_steps) 
set equ_rounds 10000 
set equ_steps  5000

# ----------------------------------------------------------- 
# save system paras 
#if 0 { 
set sys_block [open "$jobName.sys" "w"]
blockfile $sys_block write variable all
blockfile $sys_block write interactions
blockfile $sys_block write bonds all
close $sys_block 
#} 
# ----------------------------------------------------------- 
# (a) production
#if 0 {  
set equ_block [open "$jobName.equ" "w"] 
set equ_log [open "$jobName.log" "w"] 
# 
for {set i 0} {$i < $equ_rounds} {incr i} { 
    # integrate 
    integrate $equ_steps 
    # print kinetic and potential energy 
    set kin [analyze energy kinetic]  
    set tot [analyze energy total]  
    puts $equ_log [format "equ: %5d kin: %+.5e tot: %+.5e" $i $kin $tot] 
    puts [format "equ: %5d" $i] 
    flush stdout 
    # save configuration 
    blockfile $equ_block write particles "id pos type" all  
} 
close $equ_log 
close $equ_block 
#}
# ----------------------------------------------------------- 
# (b) backup status point  
#if 0 { 
set fin_block [open "$jobName.fin" "w"]
blockfile $fin_block write random
blockfile $fin_block write bitrandom
blockfile $fin_block write particles "id pos type q mol v f fix" all
close $fin_block 
sort_particles 
#} 

# ending time 
set et [clock format [clock seconds] -format %Y,%b.%d-%H:%M:%S] 
puts [format "%s ~ %s" $st $et]   

