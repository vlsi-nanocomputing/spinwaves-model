#!/bin/bash
echo "---------------  synt.sh starts ------------------"
echo "-Elimination of work directory and we creation of a new one"
rm -rf work
mkdir work
echo "-Inizialization of synopsys" 

# source /software/scripts/init_synopsys_64.18
echo "-Execution of synthesys through synt.scr, output is redirected in synt_out/dc_shell_output.txt"
dc_shell-t -f synt.scr > ./synt_out/dc_shell_output.txt
echo -e "--------------- end\n\n"
