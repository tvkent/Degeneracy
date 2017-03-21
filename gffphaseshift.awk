#!/bin/bash
#adjust feature bed file based on phase

BEGIN {OFS="\t"}
{
	if ($4=="+" && $5=="1"){print $1,($2+1),$3,".",".",$4}
	else if ($4=="+" && $5=="2"){print $1,($2+2),$3,".",".",$4}
	else if ($4=="-" && $5=="1"){print $1,$2,($3-1),".",".",$4}
	else if ($4=="-" && $5=="2"){print $1,$2,($3-2),".",".",$4}
	else {print $1,$2,$3,".",".",$4}
}

