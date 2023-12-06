#!/bin/bash
kpoints=''$IRVSPDATA'/max_KPOINTS_VASP'

if [ X$1 = X ];then
  echo "!!!!!! ./IRphx.sh  #sgn  !!!!!! "
  exit
fi

nks=`sed -n '2p' $kpoints/KPOINTS_$1.txt`
edks=`echo "$nks+3"|bc`
nk=(`sed -n '4,'$edks'p' $kpoints/KPOINTS_$1.txt| awk  '{print $1,$2,$3}'`)

if [ ! -d "./tmp" ];then
   echo "./tmp/ not found !!!"
   exit
fi
bb1=(`grep '<b1>' ./tmp/pwscf.xml |head -1|awk -F '>' '{print $2}'|awk -F '<' '{print $1}'|awk '{print $0}'`)
bb2=(`grep '<b2>' ./tmp/pwscf.xml |head -1|awk -F '>' '{print $2}'|awk -F '<' '{print $1}'|awk '{print $0}'`)
bb3=(`grep '<b3>' ./tmp/pwscf.xml |head -1|awk -F '>' '{print $2}'|awk -F '<' '{print $1}'|awk '{print $0}'`)
b1=(`printf "%10.6f%10.6f%10.6f" "${bb1[@]}"`)
b2=(`printf "%10.6f%10.6f%10.6f" "${bb2[@]}"`)
b3=(`printf "%10.6f%10.6f%10.6f" "${bb3[@]}"`)
echo $0 $1


if [ -f "ph_wf.dat" ];then
  rm ph_wf.dat
fi

ldisp=`sed -n '/ldisp.*/p' ph.inp`
if [ -n "$ldisp" ];then
  sed -i 's/ldisp.*/ldisp=.false./g' ph.inp
fi

for((i=0;i<"$nks";i++))
do

if [ -f "q$i.dyn" ];then
   echo "q$i.dyn is added to ph_wf.dat successfully !!!"
   if [ $i+1=="$nks" ]; then
      python ./dyn2wf.py  "$nks" 
   fi
else
   if [ ! -f "ph.inp" ];then
      echo "ph.inp not found !!!"
      exit
   fi
   cp ph.inp q"$i".inp 
   sed -i "s/fildyn=.*/fildyn='q$i.dyn',/g" q$i.inp
   sed -i '/^ \//,$d' q$i.inp
   printf " /\n" >> q$i.inp
   c1=`echo "${nk[3*i]}*${b1[0]}+${nk[3*i+1]}*${b2[0]}+${nk[3*i+2]}*${b3[0]}"|bc`
   c2=`echo "${nk[3*i]}*${b1[1]}+${nk[3*i+1]}*${b2[1]}+${nk[3*i+2]}*${b3[1]}"|bc`
   c3=`echo "${nk[3*i]}*${b1[2]}+${nk[3*i+1]}*${b2[2]}+${nk[3*i+2]}*${b3[2]}"|bc`
   printf "%15.8f%15.8f%15.8f\n" $c1 $c2 $c3 >> q$i.inp
   echo "q$i.inp is ready. Please run \"ph.x < q$i.inp > q$i-out\" "
fi

done
