#!/bin/sh
for f in ex[1-9].out ex1[0-2].out ex13[abc].out ex1[4-8].out \
         ex[25789].sel ex1[0245].sel ex13[abc].sel \
	 ex6A-B.sel ex6C.sel ex11adv.sel ex11trn.sel
do
   rm -f $f
done
rm -f check.out *.log *.out *.sel
