signature_list=("SBS1" "SBS2" "SBS3" "SBS4" "SBS5" "SBS6" "SBS7a" "SBS7b" "SBS7c" "SBS7d" "SBS8" "SBS9" "SBS10a" "SBS10b" "SBS10c" "SBS10d" "SBS11" "SBS12" "SBS13" "SBS14" "SBS15" "SBS16" "SBS17a" "SBS17b" "SBS18" "SBS19" "SBS20" "SBS21" "SBS22a" "SBS22b" "SBS23" "SBS24" "SBS25" "SBS26" "SBS28" "SBS29" "SBS30" "SBS31" "SBS32" "SBS33" "SBS34" "SBS35" "SBS36" "SBS37" "SBS38" "SBS39" "SBS40a" "SBS40b" "SBS40c" "SBS41" "SBS42" "SBS44" "SBS84" "SBS85" "SBS86" "SBS87" "SBS88" "SBS89" "SBS90" "SBS91" "SBS92" "SBS93" "SBS94" "SBS96" "SBS97" "SBS98" "SBS99")

parallel -j20 python3 /home/maria/cactus_target_size/scripts/3_multiply_matrices_find_target_size_tsb.py ::: "${signature_list[@]}"




