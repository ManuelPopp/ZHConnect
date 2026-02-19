seq 0 24 | xargs -n 1 -P 5 -I {} \
bash -c 'sleep $(( {} * 8 )); nohup python3 /lud11/poppman/shared/dami/py3/resultsECA.py \
-rst "top9_C_comp_gridbased.tif" \
-residx {} > /lud11/poppman/shared/dami/log/nohup{}.out 2> /lud11/poppman/shared/dami/log/nohup{}.err'

seq 0 24 | xargs -n 1 -P 4 -I {} \
bash -c 'sleep $(( {} * 4 )); nohup python3 /lud11/poppman/shared/dami/py3/resultsECA.py \
-rst "top9_C_comp_patchbased.tif" \
-residx {} > /lud11/poppman/shared/dami/log/nohup{}.out 2> /lud11/poppman/shared/dami/log/nohup{}.err'

seq 0 24 | xargs -n 1 -P 4 -I {} \
bash -c 'sleep $(( {} * 4 )); nohup python3 /lud11/poppman/shared/dami/py3/resultsECA.py \
-rst "top9_C_con_gridbased.tif" \
-residx {} > /lud11/poppman/shared/dami/log/nohup{}.out 2> /lud11/poppman/shared/dami/log/nohup{}.err'

seq 0 24 | xargs -n 1 -P 4 -I {} \
bash -c 'sleep $(( {} * 4 )); nohup python3 /lud11/poppman/shared/dami/py3/resultsECA.py \
-rst "top9_C_con_patchbased.tif" \
-residx {} > /lud11/poppman/shared/dami/log/nohup{}.out 2> /lud11/poppman/shared/dami/log/nohup{}.err'

seq 0 24 | xargs -n 1 -P 5 -I {} \
bash -c 'sleep $(( {} * 8 )); nohup python3 /lud11/poppman/shared/dami/py3/resultsECA.py \
-rst "top9_R_comp_gridbased.tif" \
-residx {} > /lud11/poppman/shared/dami/log/nohup{}.out 2> /lud11/poppman/shared/dami/log/nohup{}.err'

seq 0 24 | xargs -n 1 -P 4 -I {} \
bash -c 'sleep $(( {} * 4 )); nohup python3 /lud11/poppman/shared/dami/py3/resultsECA.py \
-rst "top9_R_comp_patchbased.tif" \
-residx {} > /lud11/poppman/shared/dami/log/nohup{}.out 2> /lud11/poppman/shared/dami/log/nohup{}.err'

seq 0 24 | xargs -n 1 -P 4 -I {} \
bash -c 'sleep $(( {} * 4 )); nohup python3 /lud11/poppman/shared/dami/py3/resultsECA.py \
-rst "top9_R_con_gridbased.tif" \
-residx {} > /lud11/poppman/shared/dami/log/nohup{}.out 2> /lud11/poppman/shared/dami/log/nohup{}.err'

seq 0 24 | xargs -n 1 -P 4 -I {} \
bash -c 'sleep $(( {} * 4 )); nohup python3 /lud11/poppman/shared/dami/py3/resultsECA.py \
-rst "top9_R_con_patchbased.tif" \
-residx {} > /lud11/poppman/shared/dami/log/nohup{}.out 2> /lud11/poppman/shared/dami/log/nohup{}.err'

seq 0 24 | xargs -n 1 -P 4 -I {} \
bash -c 'sleep $(( {} * 4 )); nohup python3 /lud11/poppman/shared/dami/py3/resultsECA.py \
-rst "baseline_gridbased.tif" \
-residx {} > /lud11/poppman/shared/dami/log/nohup{}.out 2> /lud11/poppman/shared/dami/log/nohup{}.err'
