#!/bin/bash
#to be run as "user@Mayu$ sh var/test.sh"

echo REGRESSION TEST FOR MAYU

./Mayu.pl -B example.csv -C tardecdb.fa -m testoutput -verbose

ref="protFDR
0.04038197"
tgt=$(awk 'NR==1||$3==0.01{print $20}' testoutput_main_1.08.txt)

if [ "$ref" == "$tgt" ]; then
	echo TEST OK
else
	echo TEST FAIL
fi

rm testoutput_main_*
