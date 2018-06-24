less 1node.papi.2493170.out | grep AE | sed 's/AE: //g'
less 1node.papi.2493170.out | grep t_total | cut -d ',' -f 2 | cut -c 10-
