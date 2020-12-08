#!/bin/bash


N=1
for (( i=1; i<=5; i++ ))
do
	let N=N*2

	cd ../VHDL
	while read line
	do
		if [[ "$line" =~ 'NATURAL' ]];
		then
			echo "GENERIC(N : NATURAL:=$N);" >> RCA_nbit_new.vhd
		else
			echo $line >> RCA_nbit_new.vhd
		fi
	done < RCA_nbit.vhd

	rm RCA_nbit.vhd
	mv RCA_nbit_new.vhd RCA_nbit.vhd	
	cd ../synt


	while read line
	do
		if [[ "$line" =~ 'ext' ]];
		then
			echo "variable ext _${N}bit" >> variables_new.scr
		else
			echo $line >> variables_new.scr
		fi
	done < variables.scr

	rm variables.scr
	mv variables_new.scr variables.scr	

	

	echo "----------------- ${N}-bit RCA analysis -------------"
source synt_save.sh
	
done
