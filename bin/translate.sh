faa=$(echo $1 | sed -E "s/(.*)\..*/\1/").faa
fastatranslate $1 -F 1 | sed "/>/ s/ \\[translate(1)\\]//" > $faa.tmp
# linearize
awk '
	# credits to https://gist.github.com/lindenb/2c0d4e11fd8a96d4c345
	/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;}
	     {printf("%s",$0);}
	END  {printf("\n");}
' $faa.tmp > $faa
rm $faa.tmp
