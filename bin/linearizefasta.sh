#!/bin/bash
awk '
	# credits to https://gist.github.com/lindenb/2c0d4e11fd8a96d4c345
	/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;}
	     {printf("%s",$0);}
	END  {printf("\n");}
' $1
