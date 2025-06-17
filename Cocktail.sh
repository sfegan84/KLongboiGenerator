#Coctail generator
#number of hours of running to simulate
Hr=$1;

rm -rf ./tmp
mkdir tmp

htmp=$((Hr*144082));

echo $Hr " " $htmp 

KLGenerator_hddm_V3 -M${htmp} -Fgenerated.root -Ekaon:histo:0.0:10.0 -Rkl1 -SsolC

mv generated.hddm ./tmp/kl1.hddm
mv generated.root ./tmp/kl1.root

#kl2
htmp=$((Hr*8141));

echo $Hr " " $htmp 

KLGenerator_hddm_V3 -M${htmp} -Fgenerated.root -Ekaon:histo:0.0:10.0 -Rkl2 -SsolC

mv generated.hddm ./tmp/kl2.hddm
mv generated.root ./tmp/kl2.root

#kl3
htmp=$((Hr*235));

echo $Hr " " $htmp 

KLGenerator_hddm_V3 -M${htmp} -Fgenerated.root -Ekaon:histo:0.0:10.0 -Rkl3 -SsolC

mv generated.hddm ./tmp/kl3.hddm
mv generated.root ./tmp/kl3.root

#kl4
htmp=$((Hr*14479));

echo $Hr " " $htmp 

KLGenerator_hddm_V3 -M${htmp} -Fgenerated.root -Ekaon:histo:0.0:10.0 -Rkl4 -SsolC

mv generated.hddm ./tmp/kl4.hddm
mv generated.root ./tmp/kl4.root

#kl5
htmp=$((Hr*6905));

echo $Hr " " $htmp 

KLGenerator_hddm_V3 -M${htmp} -Fgenerated.root -Ekaon:histo:0.0:10.0 -Rkl5 -SsolC

mv generated.hddm ./tmp/kl5.hddm
mv generated.root ./tmp/kl5.root

#kl6
htmp=$((Hr*6905));

echo $Hr " " $htmp 

KLGenerator_hddm_V3 -M${htmp} -Fgenerated.root -Ekaon:histo:0.0:10.0 -Rkl6 -SsolC

mv generated.hddm ./tmp/kl6.hddm
mv generated.root ./tmp/kl6.root

#kl9
htmp=$((Hr*152224));

echo $Hr " " $htmp 

KLGenerator_hddm_V3 -M${htmp} -Fgenerated.root -Ekaon:histo:0.0:10.0 -Rkl9 -SsolC

mv generated.hddm ./tmp/kl9.hddm
mv generated.root ./tmp/kl9.root

rm ./Cocktail_kl.root
hadd ./Cocktail_kl.root ./tmp/kl1.root ./tmp/kl2.root ./tmp/kl3.root ./tmp/kl4.root ./tmp/kl5.root ./tmp/kl6.root  ./tmp/kl9.root 

hddm_merge_files -oCocktail_kl.hddm ./tmp/kl1.hddm  ./tmp/kl2.hddm  ./tmp/kl3.hddm  ./tmp/kl4.hddm  ./tmp/kl5.hddm  ./tmp/kl6.hddm  ./tmp/kl9.hddm
#done

#num1=10
#num2=20

#ans=$((num1 * num2))

#echo $ans
