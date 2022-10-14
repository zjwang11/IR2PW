import numpy as np
from numpy.linalg import *
import os
from datetime import datetime


def mode():
    if os.path.exists('./ph.fc') is not True or os.path.exists('./matdyn.inp') is not True:
        print("We need ph.fc & matdyn.inp...")
    else:
        if os.path.exists('./matdyn.modes') is not True:
            exist()
        if os.path.exists('./wfc') == 0 :
            os.mkdir('./wfc')
        
        f1 = open('./ph.fc')
        fcdata = f1.readline().split()
        global fcnat
        fcnat = int(fcdata[1])
        f1.close()

        f2 = open('./matdyn.inp')
        lines = f2.readlines()
        for i in range(len(lines)):
            if '/' in lines[i]:
                kps = int(lines[i+1])
                kpn = i+1
        nq = 0
        for i in range(kps-1):
            nq += int(lines[kpn+i+1].split()[3])
        nq += 1
        f2.close()


        eqmode = 3*fcnat*(fcnat+1)+5             # modelines = eqmode * nq  # num of mode lines

        u = np.zeros([nq,3*fcnat,6*fcnat])
        fr = np.zeros([nq,3*fcnat,1])
        f3 = open('./matdyn.modes')
        mlines = f3.readlines()
        f3.close()
        
        pfile = open('./wfc/modes.dat','w')
        for i in range(nq):
            for j in range(3*fcnat):
                fr[i][j][0] = float(mlines[eqmode*i+4+j*(fcnat+1)].split("=")[2].split()[0])
                for k in range(fcnat):
                    for l in range(6):
                            u[i][j][k*6+l] = float(mlines[eqmode*i+5+j*(fcnat+1)+k].split()[l+1])
       
            
        print('{:8d}{:8d}{:8d}'.format(nq,3*fcnat,3*fcnat),file=pfile)
        for i in range(nq):
            for j in range(3*fcnat):
                print('{:13.6f}'.format(fr[i][j][0]),file=pfile)
                
                for kk in range((3*fcnat)//5):
                    xxx = (kk+1)*5
                    for kkk in range(xxx-5, xxx):
                        print("{:12.6f}{:12.6f}".format(u[i][j][2*kkk],u[i][j][2*kkk+1]), end='',file=pfile)
                    print(' ',file=pfile)
                if 3*fcnat % 5 > 0:
                    for k in range(xxx, xxx+(3*fcnat)%5):
                        print("{:12.6f}{:12.6f}".format(u[i][j][2*k],u[i][j][2*k+1]), end='',file=pfile)
                    print(' ',file=pfile)
        pfile.close()
        f4 = open('./wfc/lda_hr.dat','w')
        line=" Generated on "+str(datetime.now())+"\n"
        f4.write(line)
        f4.write('{:8d}\n'.format(3*fcnat))
        f4.write('{:8d}\n'.format(1))
        f4.write('{:8d}\n'.format(1))
        for i in range(1,3*fcnat+1):
            for j in range(1,3*fcnat+1):
                f4.write('{:8d}{:8d}{:8d}{:8d}{:8d}{:20.10f}{:20.10f}\n'.format(0,0,0,j,i,0,0))
        f4.close()
#        exist()

        
        

#def bcskpoints():
#    f0 = open('./ph.fc')
#    l = f0.readline().split()
#    ntyp = int(l[0])
#    natm = int(l[1])
#    brav = int(l[2])
#    cell_para = np.zeros([3,3])
#    atom_name = np.zeros([ntyp],dtype=object)
#    label_at = np.zeros([ntyp],dtype=int)
#    atoms = np.zeros([ntyp],dtype=int)
#    atom_posi = np.zeros([natm,3])
#    if brav == 0 :
#        for i in range(3):
#            l = f0.readline().split()
#            for j in range(3):
#                cell_para[i][j] = l[j]
#    for i in range(ntyp):
#        l = f0.readline().split('\'')
#        atom_name[i] = l[1]
#        label_at[i]  = int(l[0])
#    for i in range(natm):
#        l = f0.readline().split()
#        for j in range(ntyp):
#            if int(l[1]) == label_at[j]:
#                atoms[j] += 1
#        for k in range(3):
#            atom_posi[i][k] = l[k+2]
#    f0.close()
#    f = open('./fc.vasp','w')
#    f.write('ph.fc write \n')
#    f.write('1.0\n')
#    for i in range(3):
#        f.write('{:18.10f}{:18.10f}{:18.10f}\n'.format(cell_para[i][0],cell_para[i][1],cell_para[i][2]))
#    for i in range(ntyp):
#        f.write('{:5s}'.format(atom_name[i]))
#    f.write('\n')
#    for i in range(ntyp):
#        f.write('{:<5d}'.format(atoms[i]))
#    f.write('\n')
#    f.write('Cartesian\n')
#    for i in range(natm):
#        f.write('{:18.10f}{:18.10f}{:18.10f}\n'.format(atom_posi[i][0],atom_posi[i][1],atom_posi[i][2]))
#    f.write('\n')
#    f.close()
#    os.system(r'phonopy --symmetry --tolerance 0.01 -c fc.vasp > symmetry')
#    global sgn
#    sgn = os.popen('grep "space_group_number:" symmetry|awk \'{print $2}\'|xargs echo -n').read()
#    
#    os.system('rm BPOSCAR PPOSCAR symmetry fc.vasp')
#    address = '/home/Zhangrh/max_KPOINTS_VASP/'
#    f1 = open(address+'KPOINTS_'+sgn+'.txt')
#    kp=f1.readlines()
#    nkps = int(kp[1])
#    rec = np.zeros([nkps,3])
#    for i in range(nkps):
#        kpsf = kp[3+i].split()
#        for j in range(3):
#            rec[i][j]=float(kpsf[j])
#    f1.close()
#    f2 = open('./qeinput_{}'.format(sgn),'w')
#    print('/',file=f2)
#    print(nkps,file=f2)
#    for i in range(nkps):
#        print("{:14.8f}{:14.8f}{:14.8f}{:8d}".format(rec[i][0],rec[i][1],rec[i][2],1),file=f2)
#    f2.close()
#     
#
#def exist():
#    bcskpoints()
#    os.system("sed -i '/\//,$d' matdyn.inp")
#    os.system("cat qeinput_"+sgn+" >> matdyn.inp")
#    os.system('matdyn.x -in matdyn.inp > matdyn.out')
     


mode()

