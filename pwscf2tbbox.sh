#!/bin/sh
name='tbbox.in'
infile='scf.out'
kpoints=''$IRVSPDATA'/max_KPOINTS_VASP'

if [ X$1 = X ];then
  echo "!!!!!! ./pwscf2tbbox.sh  #sgn  !!!!!! "
  exit
fi

nks=`sed -n '2p' $kpoints/KPOINTS_$1.txt`
edks=`echo "$nks+3"|bc`
nk=(`sed -n '4,'$edks'p' $kpoints/KPOINTS_$1.txt| awk  '{print $1,$2,$3}'`)

if [ ! -f "$infile" ];then
   echo "$infile not found !!!"
   exit
fi
nat=`grep 'number of atoms/cell' $infile |awk '{print $5}'`
po=(`grep -A "$nat" 'positions (cryst. coord.)'  $infile |grep 'tau'| awk  '{print  $7 ,$8, $9 }'`)
atme=(`grep -A "$nat" 'positions (cryst. coord.)' $infile |grep 'tau'| awk  '{print  $2 }'`)

count=($(seq $nat))
for ((i=0;i<"$nat";i++))
do
  for ((j=i+1;j<"$nat";j++))
  do
    if [ ${atme[i]} = ${atme[j]} ];then
       count[j]="${count[i]}"
    fi
  done
done
printf "  case = ph\n" > $name
printf "\n" >> $name
printf "  proj:\n" >> $name
printf "  orbt = 1 \n" >> $name
printf "  ntau = $nat \n" >> $name
for ((i=0;i<"$nat";i++))
do
printf "%12.6f%12.6f%12.6f%5.0f 3\n" ${po[3*i]} ${po[3*i+1]} ${po[3*i+2]} ${count[i]} >> $name
done
printf "  end projections\n" >> $name
printf "\n" >> $name
printf "  kpoint:\n" >> $name
printf "  kmesh = 1 \n" >> $name
printf "  Nk = %-3d\n" $nks >> $name
for((i=0;i<"$nks";i++))
do
  printf "%15.8f%15.8f%15.8f\n" ${nk[3*i]} ${nk[3*i+1]} ${nk[3*i+2]} >> $name
done
printf "  end kpoint_path\n" >> $name
printf "\n" >> $name
printf "\n" >> $name
printf "  unit_cell:\n" >> $name
cryax=(`grep -A 3 "crystal axes:" $infile |tail -3|awk '{printf "%11s%11s%11s\n", $4 ,$5 ,$6}'`)
repax=(`grep -A 3 "reciprocal axes:" $infile |tail -3|awk '{printf "%11s%11s%11s\n", $4 ,$5 ,$6}'`)
for((i=0;i<3;i++))
do 
printf "  %11.6f%11.6f%11.6f  %11.6f%11.6f%11.6f\n" ${cryax[3*i]} ${cryax[3*i+1]} ${cryax[3*i+2]} ${repax[3*i]} ${repax[3*i+1]} ${repax[3*i+2]} >> $name
done


sym=`grep "Sym. Ops." $infile |awk '{print $1}'`
crystax=(`grep -A 3 "crystal axes:" "$infile" |tail -3|awk '{printf "%11s%11s%11s\n", $4 ,$5 ,$6}'`)

if grep -q 'Time Reversal' "$infile"; then 
  trev=(`grep 'Time Reversal' "$infile"|awk '{print $3}'`)
fi
ftau1=(`grep  'cryst.   s' "$infile"|awk -F '(' '{print $4}'|awk -F ')' '{print $1}'`)
itau=(`grep 'cryst.   s' "$infile" |grep 'f'|awk -F '(' '{print $2}'|awk -F ')' '{print $1}'`)
for i in $(seq 1 $sym)
do
  tau1[i]=0
done

i=0
for j in ${itau[@]};
do
  tau1[j-1]=`echo "-1*${ftau1[i]}"|bc`
  i=`echo "$i+1"|bc`
done


ftau2=(`sed -n '/cryst.   s/{n;p;}' "$infile"|awk -F '(' '{print $3}'|awk -F ')' '{print $1}'`)

for i in $(seq 1 $sym)
do
  tau2[i]=0
done

i=0
for j in ${itau[@]};
do
  tau2[j-1]=`echo "-1*${ftau2[i]}"|bc`
  i=`echo "$i+1"|bc`
done


ftau3=(`sed -n '/cryst.   s/{n;n;p;}' "$infile"|awk -F '(' '{print $3}'|awk -F ')' '{print $1}'`)

for i in $(seq 1 $sym)
do
  tau3[i]=0
done

i=0
for j in ${itau[@]};
do
  tau3[j-1]=`echo "-1*${ftau3[i]}"|bc`
  i=`echo "$i+1"|bc`
done


det=(`grep 'isym =' "$infile" |awk '{print $4}'`)
alph=(`grep 'isym =' "$infile" |awk '{print $4}'`)
alphi=(`grep 'isym =..*inv. ' "$infile" |awk '{print $5}'`)
typx=(`grep 'isym =' "$infile" |awk '{print $8}'`)
typxi=(`grep 'isym =' "$infile" |awk '{print $9}'`)
n_xy=(`grep 'isym =' "$infile"|awk '{print $10}'`)
n_xyi=(`grep 'isym =' "$infile"|awk '{print $11}'`)

k=0
j=0
kk=0
jj=0
for ((i=0;i<"$sym";i++))
do
if [[ ${det[i]} == inv* ]]; then
  deta[i]=-1
  if [[ ${det[i]} == inversion ]]; then
    alpha[i]=0
    typex[i]='cart.'
    n_xyz[i]='[1,0,0]'
  else
    alpha[i]=`echo "-1*${alphi[j]}"|bc` 
    typex[i]="${typxi[kk]}"
    n_xyz[i]="${n_xyi[j]}"
    j=`echo "$j+1"|bc`
    kk=`echo "$kk+1"|bc`
  fi
else
  deta[i]=1
  if [[ ${det[i]} == ident* ]]; then
    alpha[i]=0 
    typex[i]='cart.'
    n_xyz[i]='[1,0,0]'
    k=`echo "$k+1"|bc`
  else
    alpha[i]=`echo "-1*${alph[k]}"|bc` 
    typex[i]="${typx[kk]}"
    n_xyz[i]="${n_xy[jj]}"
    k=`echo "$k+1"|bc`
    kk=`echo "$kk+1"|bc`
    jj=`echo "$jj+1"|bc`
  fi
fi
done

nxyz=(`echo "${n_xyz[@]}"|sed "s/\]/ /g"|sed "s/\[/ /g" |sed "s/,/ /g"|awk '{print $0}'`)

x=0
for ((i=0;i<"$sym";i++))
do
nx[i]="${nxyz[x]}"
ny[i]="${nxyz[x+1]}"
nz[i]="${nxyz[x+2]}"
x=`echo "$x+3"|bc`
if [[ ${typex[i]} == cart. ]]; then
  n_x[i]="${nx[i]}"
  n_y[i]="${ny[i]}"
  n_z[i]="${nz[i]}"
else
  xs=`echo "(${nx[i]}*${crystax[0]}+${ny[i]}*${crystax[3]}+${nz[i]}*${crystax[6]})"|bc`
  ys=`echo "(${nx[i]}*${crystax[1]}+${ny[i]}*${crystax[4]}+${nz[i]}*${crystax[7]})"|bc`
  zs=`echo "(${nx[i]}*${crystax[2]}+${ny[i]}*${crystax[5]}+${nz[i]}*${crystax[8]})"|bc`
  nf=`echo "sqrt(($xs)^2+($ys)^2+($zs)^2)"|bc`
  n_x[i]=`echo "scale=6;$xs/$nf"|bc`
  n_y[i]=`echo "scale=6;$ys/$nf"|bc`
  n_z[i]=`echo "scale=6;$zs/$nf"|bc`
fi
  ii=`echo "$i+1"|bc`
  printf "%5s%13.6f%13.6f%13.6f%13.6f%13.6f%13.6f%13.6f%13.6f\n" $ii ${deta[i]} ${alpha[i]} ${n_x[i]} ${n_y[i]} ${n_z[i]} ${tau1[i]} ${tau2[i]} ${tau3[i]}  >> $name
done

printf "  end unit_cell_cart\n" >> $name
printf "\n" >> $name


