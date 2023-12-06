from datetime import datetime
import numpy as np
from numpy import pi, sqrt
import os, sys


######README#####
# Prepare the file "ph.fc" without NAC, if your ph.fc contain dielectric constant and Born effective charges, fc2hr.py will collect them in allborn and break. You need to remove the dielectric constant and Born effective charges data from Gamma point dynamical matrix file and do q2r.x again to get the ph.fc without NAC.

# $ python fc2hr.py $asr(optional)

# If $asr = 1, the same as QE/matdyn.x's 'simple' tag : previous implementation of the asr used(3 translational asr imposed by correction of the diagonal elements of the force constants matrix.
# Otherwise, you can obtain phhr_cm1.dat directly without considering asr.

######README#####



################  CONSTANTS ################
AMU_SI = 1.6605402e-27  #[kg]
ELECTRONMASS_SI  = 9.10938215e-31 
ELECTRONVOLT_SI  = 1.60217733e-19 #[J]
AMU_AU = AMU_SI / ELECTRONMASS_SI
AMU_RY = AMU_AU / 2.0   #amu_ry = 911.4442431086565
PlanckConstant = 4.13566733e-15   #[eV s]
Hbar = PlanckConstant / (2 * pi)  # [eV s]
H_PLANCK_EV = PlanckConstant * ELECTRONVOLT_SI
tpi = 2.0 * pi
fpi = 4.0 * pi
C_SI = 299792458  # [m/s]
Mu0 = 4.0e-7 * pi  # [Hartree/m]
Epsilon0 = 1.0 / Mu0 / C_SI**2  # [C^2 / N m^2]
HARTREE_EV = ELECTRONMASS_SI * ELECTRONVOLT_SI / 16 / pi**2 / Epsilon0**2 / Hbar**2  # Hartree [eV] 27.211398 
RYDBERG_EV = HARTREE_EV/2.0
AU_SEC = H_PLANCK_EV/tpi/HARTREE_EV
AU_PS = AU_SEC * 1.0e+12
AU_TERAHERTZ = AU_PS
RY_TO_THZ = 1.0 / AU_TERAHERTZ / fpi
RY_TO_CMM1 =1.0e+10 * RY_TO_THZ / C_SI
Angstrom = 1.0e-10  # [m]
BOHR_ANG =  4e10 * pi * Epsilon0 * Hbar**2 / ELECTRONMASS_SI  # Bohr radius [A] 0.5291772
BOHR_RADIUS_SI = BOHR_ANG * 1e-10
VaspToTHz = sqrt(ELECTRONVOLT_SI / AMU_SI) / Angstrom / (2 * pi) / 1e12  # [THz] 15.633302
THzToCm = 1.0e12 / (C_SI * 100)
zi = complex(0,1)
e2 = 2.0 
eps = 1.0e-6
################  CONSTANTS ################

################  ATOM DATA ################
atom_data = [
    [0, "X", "X", None],  # 0
    [1, "H", "Hydrogen", 1.00794],  # 1
    [2, "He", "Helium", 4.002602],  # 2
    [3, "Li", "Lithium", 6.941],  # 3
    [4, "Be", "Beryllium", 9.012182],  # 4
    [5, "B", "Boron", 10.811],  # 5
    [6, "C", "Carbon", 12.0107],  # 6
    [7, "N", "Nitrogen", 14.0067],  # 7
    [8, "O", "Oxygen", 15.9994],  # 8
    [9, "F", "Fluorine", 18.9984032],  # 9
    [10, "Ne", "Neon", 20.1797],  # 10
    [11, "Na", "Sodium", 22.98976928],  # 11
    [12, "Mg", "Magnesium", 24.3050],  # 12
    [13, "Al", "Aluminium", 26.9815386],  # 13
    [14, "Si", "Silicon", 28.0855],  # 14
    [15, "P", "Phosphorus", 30.973762],  # 15
    [16, "S", "Sulfur", 32.065],  # 16
    [17, "Cl", "Chlorine", 35.453],  # 17
    [18, "Ar", "Argon", 39.948],  # 18
    [19, "K", "Potassium", 39.0983],  # 19
    [20, "Ca", "Calcium", 40.078],  # 20
    [21, "Sc", "Scandium", 44.955912],  # 21
    [22, "Ti", "Titanium", 47.867],  # 22
    [23, "V", "Vanadium", 50.9415],  # 23
    [24, "Cr", "Chromium", 51.9961],  # 24
    [25, "Mn", "Manganese", 54.938045],  # 25
    [26, "Fe", "Iron", 55.845],  # 26
    [27, "Co", "Cobalt", 58.933195],  # 27
    [28, "Ni", "Nickel", 58.6934],  # 28
    [29, "Cu", "Copper", 63.546],  # 29
    [30, "Zn", "Zinc", 65.38],  # 30
    [31, "Ga", "Gallium", 69.723],  # 31
    [32, "Ge", "Germanium", 72.64],  # 32
    [33, "As", "Arsenic", 74.92160],  # 33
    [34, "Se", "Selenium", 78.96],  # 34
    [35, "Br", "Bromine", 79.904],  # 35
    [36, "Kr", "Krypton", 83.798],  # 36
    [37, "Rb", "Rubidium", 85.4678],  # 37
    [38, "Sr", "Strontium", 87.62],  # 38
    [39, "Y", "Yttrium", 88.90585],  # 39
    [40, "Zr", "Zirconium", 91.224],  # 40
    [41, "Nb", "Niobium", 92.90638],  # 41
    [42, "Mo", "Molybdenum", 95.96],  # 42
    [43, "Tc", "Technetium", None],  # 43
    [44, "Ru", "Ruthenium", 101.07],  # 44
    [45, "Rh", "Rhodium", 102.90550],  # 45
    [46, "Pd", "Palladium", 106.42],  # 46
    [47, "Ag", "Silver", 107.8682],  # 47
    [48, "Cd", "Cadmium", 112.411],  # 48
    [49, "In", "Indium", 114.818],  # 49
    [50, "Sn", "Tin", 118.710],  # 50
    [51, "Sb", "Antimony", 121.760],  # 51
    [52, "Te", "Tellurium", 127.60],  # 52
    [53, "I", "Iodine", 126.90447],  # 53
    [54, "Xe", "Xenon", 131.293],  # 54
    [55, "Cs", "Caesium", 132.9054519],  # 55
    [56, "Ba", "Barium", 137.327],  # 56
    [57, "La", "Lanthanum", 138.90547],  # 57
    [58, "Ce", "Cerium", 140.116],  # 58
    [59, "Pr", "Praseodymium", 140.90765],  # 59
    [60, "Nd", "Neodymium", 144.242],  # 60
    [61, "Pm", "Promethium", None],  # 61
    [62, "Sm", "Samarium", 150.36],  # 62
    [63, "Eu", "Europium", 151.964],  # 63
    [64, "Gd", "Gadolinium", 157.25],  # 64
    [65, "Tb", "Terbium", 158.92535],  # 65
    [66, "Dy", "Dysprosium", 162.500],  # 66
    [67, "Ho", "Holmium", 164.93032],  # 67
    [68, "Er", "Erbium", 167.259],  # 68
    [69, "Tm", "Thulium", 168.93421],  # 69
    [70, "Yb", "Ytterbium", 173.054],  # 70
    [71, "Lu", "Lutetium", 174.9668],  # 71
    [72, "Hf", "Hafnium", 178.49],  # 72
    [73, "Ta", "Tantalum", 180.94788],  # 73
    [74, "W", "Tungsten", 183.84],  # 74
    [75, "Re", "Rhenium", 186.207],  # 75
    [76, "Os", "Osmium", 190.23],  # 76
    [77, "Ir", "Iridium", 192.217],  # 77
    [78, "Pt", "Platinum", 195.084],  # 78
    [79, "Au", "Gold", 196.966569],  # 79
    [80, "Hg", "Mercury", 200.59],  # 80
    [81, "Tl", "Thallium", 204.3833],  # 81
    [82, "Pb", "Lead", 207.2],  # 82
    [83, "Bi", "Bismuth", 208.98040],  # 83
    [84, "Po", "Polonium", None],  # 84
    [85, "At", "Astatine", None],  # 85
    [86, "Rn", "Radon", None],  # 86
    [87, "Fr", "Francium", None],  # 87
    [88, "Ra", "Radium", None],  # 88
    [89, "Ac", "Actinium", None],  # 89
    [90, "Th", "Thorium", 232.03806],  # 90
    [91, "Pa", "Protactinium", 231.03588],  # 91
    [92, "U", "Uranium", 238.02891],  # 92
    [93, "Np", "Neptunium", None],  # 93
    [94, "Pu", "Plutonium", None],  # 94
    [95, "Am", "Americium", None],  # 95
    [96, "Cm", "Curium", None],  # 96
    [97, "Bk", "Berkelium", None],  # 97
    [98, "Cf", "Californium", None],  # 98
    [99, "Es", "Einsteinium", None],  # 99
    [100, "Fm", "Fermium", None],  # 100
    [101, "Md", "Mendelevium", None],  # 101
    [102, "No", "Nobelium", None],  # 102
    [103, "Lr", "Lawrencium", None],  # 103
    [104, "Rf", "Rutherfordium", None],  # 104
    [105, "Db", "Dubnium", None],  # 105
    [106, "Sg", "Seaborgium", None],  # 106
    [107, "Bh", "Bohrium", None],  # 107
    [108, "Hs", "Hassium", None],  # 108
    [109, "Mt", "Meitnerium", None],  # 109
    [110, "Ds", "Darmstadtium", None],  # 110
    [111, "Rg", "Roentgenium", None],  # 111
    [112, "Cn", "Copernicium", None],  # 112
    [113, "Uut", "Ununtrium", None],  # 113
    [114, "Uuq", "Ununquadium", None],  # 114
    [115, "Uup", "Ununpentium", None],  # 115
    [116, "Uuh", "Ununhexium", None],  # 116
    [117, "Uus", "Ununseptium", None],  # 117
    [118, "Uuo", "Ununoctium", None],  # 118
]
################  ATOM DATA ################

def set_simple_asr(hop,nat,ntot):
    for i in range(3):
        for j in range(3):
            for na in range(nat):
                sumh = 0
                for nb in range(nat):
                    for n in range(ntot):
                        sumh += hop[n+(ntot*nb)+(ntot*nat*na)+(ntot*nat*nat*j)+(ntot*nat*nat*3*i)]
                hop[0+(ntot*na)+(ntot*nat*na)+(ntot*nat*nat*j)+(ntot*nat*nat*3*i)] -= sumh          
    return hop

def set_permutation_symmetry(hop,nat,ntot):
    hop_ = hop.copy()
    hop_hermite = hop.copy()
    for x1 in range(3):
        for x2 in range(3):
            for x1_ in range(3):
                for x2_ in range(3):
                    if x1 == x2_ and x2 == x1_ :
                        for i in range(nat):
                            for j in range(nat):
                                for i_ in range(nat):
                                    for j_ in range(nat):
                                        if i == j_ and j == i_ :
                                            ori = ntot*(j+nat*(i+nat*(x2+3*x1)))
                                            orit = ntot*(i+nat*(j+nat*(x1+3*x2)))
                                            for sc in range(ntot):
                                                hop_hermite[sc+ori] = (hop[sc+ori]+hop_[sc+orit])/2.0
    return hop_hermite

################  FORMAT ################
# for n in atom2:
#     for m in atom1:
#         < x y z  m n Re Im >  <mR|H|n0>
### D_{aa'}(k)=\sum{Rt'-Rt}{frac{phi_{t,a;t',a'}}{\sqrt{m_a*m_a'}}*exp{j*[R_t'-R_t]*k}
################  FORMAT ################


################  read fc ################
if os.path.exists('./ph.fc') is not True :
    print("ERROR: We need ph.fc !!!")
    exit(0)
cm1 = (1.E+10/C_SI)**2*(RYDBERG_EV/(BOHR_ANG)**2)
#thz = cm1 * (1/THzToCm)**2 
scale = np.array([[cm1],['cm1']],dtype=object)
k = open('./ph.fc')
fi = k.readline().split()
ntyp = int(fi[0])
nat = int(fi[1])
ibrav = fi[2]
celldm = np.zeros(6)
for i in range(6):
    celldm[i] = fi[3+i]
A = np.zeros([3,3])
for i in range(3):
    ll = k.readline().split()
    for j in range(3):
        A[i][j] = ll[j]
num = np.zeros([ntyp,2],dtype=object) 
for i in range(ntyp):
    fi = k.readline().split()
    num[i,0] = fi[0]
    num[i,1] = fi[1].split("'")[1]
for n in range(len(atom_data)):
    for i in range(ntyp):
        if atom_data[n][1] == num[i,1]:
            num[i,1] = atom_data[n][3]

op = np.ones([3,3]) ###xyz orb
atom_mass = np.zeros([nat,nat])
num_mass = np.zeros([1,nat])   
cart_pos = np.zeros([nat,3])
cry_pos = np.zeros([nat,3])
for i in range(nat):
    fi = k.readline().split()
    mat_vec = np.mat([[A[0][0],A[1][0],A[2][0]],[A[0][1],A[1][1],A[2][1]],[A[0][2],A[1][2],A[2][2]]],float)
    mat_pos = np.mat([fi[2],fi[3],fi[4]],float).T
    res = np.linalg.solve(mat_vec,mat_pos)
    for l in range(3):
        cart_pos[i][l] = fi[l+2]
        cry_pos[i][l] = res[l]
    for j in range(ntyp):
        if float(fi[1]) == float(num[j][0]):
            num_mass[0][i] = num[j,1]
for i in range(nat):
    for j in range(nat):
        atom_mass[i,j] = np.sqrt(num_mass[0,i]*num_mass[0,j])
atom_obmass = np.kron(atom_mass,op)
at_obmass = 1/atom_obmass
fi = k.readline().split()

####### write epsilon and born effective charge  #####
if fi[0] == 'T':
    print('Please remove the dielectric constant and Born effective charges data from the Gamma point dyn file and run q2r.x again to get correct ph.fc!!!')
    born=open('./allborn','w')
    epsilon=np.zeros([3,3])
    zstar=np.zeros([nat,3,3])
    for i in range(3):
        epsilon[i,:] = k.readline().split()
        born.write("{:24.12f}{:24.12f}{:24.12f}\n".format(epsilon[i,0],epsilon[i,1],epsilon[i,2],))
    for i in range(nat):
        y = k.readline()
        zstar[i,0,:]=k.readline().split()
        zstar[i,1,:]=k.readline().split()
        zstar[i,2,:]=k.readline().split()
        born.write('  {}  #{}\n'.format('atom',i+1))
        born.write("{:15.7f}{:15.7f}{:15.7f}\n".format(zstar[i,0,0],zstar[i,0,1],zstar[i,0,2],))
        born.write("{:15.7f}{:15.7f}{:15.7f}\n".format(zstar[i,1,0],zstar[i,1,1],zstar[i,1,2],))
        born.write("{:15.7f}{:15.7f}{:15.7f}\n".format(zstar[i,2,0],zstar[i,2,1],zstar[i,2,2],))
    born.close()
    exit()
####### end epsilon and born effective charge  #####

y = k.readline()
nq = y.split()
nks=[]
for i in range(0,3):
    nks.append(nq[i])
ntot=int(nks[0]) * int(nks[1]) * int(nks[2])
x = k.read()
k.close()
list=x.split()
a=[]  ##atom1
ax=[] ##polar1xyz
b=[]  ##atom2
bx=[] ##polar2xyz
titl=(len(list))/4  ##line_number
nat_l = int(titl/(ntot+1)) ##each_atoms_xyzs_count  {(x,y,z)1} {(x,y,z)2} {sum(atom1)} {sum(atom2)}
bas = int(np.sqrt(nat_l))
for i in range(0,nat_l):
    ax.append(int(list[0+(4*(ntot+1)*i)]))
    bx.append(int(list[1+(4*(ntot+1)*i)]))
    a.append( int(list[2+(4*(ntot+1)*i)]))
    b.append( int(list[3+(4*(ntot+1)*i)]))
################  end read fc ################

################  find multi ################
nr = np.zeros(3)
for i in range(3):
    nr[i] = float(nks[i])
rx=int(nr[0])
ry=int(nr[1])
rz=int(nr[2])
sp_pos = np.zeros([nat*ntot,3])
for i in range(nat): 
    for k3 in range(rz):
        for k2 in range(ry):
            for k1 in range(rx):
                sp_pos[k1+rx*(k2+ry*(k3+i*rz))] = cry_pos[i] + [k1,k2,k3]
                sp_dis = sum([c*c for c in [k1,k2,k3]])
                xsp = [k1-rx,k2,k3]
                xsp_dis = sum([c*c for c in xsp])
                if xsp_dis < sp_dis:
                    sp_pos[k1+rx*(k2+ry*(k3+i*rz))] = sp_pos[k1+rx*(k2+ry*(k3+i*rz))] - [k1,k2,k3] + xsp
                ysp = [k1,k2-ry,k3] 
                ysp_dis = sum([c*c for c in ysp])
                if ysp_dis < sp_dis:
                    sp_pos[k1+rx*(k2+ry*(k3+i*rz))] = sp_pos[k1+rx*(k2+ry*(k3+i*rz))] - [k1,k2,k3] + ysp
                zsp = [k1,k2,k3-rz]
                zsp_dis = sum([c*c for c in zsp])
                if zsp_dis < sp_dis:
                    sp_pos[k1+rx*(k2+ry*(k3+i*rz))] = sp_pos[k1+rx*(k2+ry*(k3+i*rz))] - [k1,k2,k3] + zsp
sp_b = np.zeros([nat*nat*ntot,3])
for i in range(nat):
    for j in range(nat*ntot):
        for k in range(3):
            sp_b[i*nat*ntot+j][k] = sp_pos[j][k] - cry_pos[i][k]
sp_bc = np.dot(sp_b,A)
super_pos = np.zeros([nat*ntot*27,3])
for z in range(3):
    for y in range(3):
        for x in range(3):
            for i in range(nat*ntot):
                super_pos[i+nat*ntot*(x+3*(y+3*z))] = sp_pos[i] + [rx*(x-1),ry*(y-1),rz*(z-1)]
super_b = np.zeros([27*nat*nat*ntot,3])
for a in range(27):
    for i in range(nat):
        for j in range(nat*ntot):
            for k in range(3):
                super_b[nat*ntot*(i+nat*a)+j][k] = super_pos[a*nat*ntot+j][k] - cry_pos[i][k]
super_bc = np.dot(super_b,A)
sp_posb = np.zeros([nat*nat*ntot,3])
super_posb = np.zeros([27*nat*nat*ntot,3])
for i in range(nat):
    for j in range(nat*ntot):
        sp_posb[j+i*nat*ntot] = sp_pos[j] - cry_pos[j//(ntot)]
for a in range(27):
    for i in range(nat):
        for j in range(nat*ntot):
            super_posb[nat*ntot*(nat*a+i)+j] = super_pos[a*nat*ntot+j] -cry_pos[j//(ntot)]
multi = np.zeros([nat*ntot,nat],dtype=int)
for i in range(nat*nat*ntot):
    sp_dis = sum([c*c for c in sp_bc[i]])
    for j in range(27):
        super_dis = sum([c*c for c in super_bc[i+j*nat*nat*ntot]])
        if super_dis - sp_dis < 0.0 :
            sp_dis = super_dis
            sp_bc[i] = super_bc[i+j*nat*nat*ntot]
            sp_posb[i] = super_posb[i+j*nat*nat*ntot]
    for k in range(27):
        super_dis = sum([c*c for c in super_bc[i+k*nat*nat*ntot]])        
        if abs(super_dis - sp_dis) < 0.00001:
            multi[i%(nat*ntot)][i//(nat*ntot)] += 1
################  end find multi ################

################  change order ################
########  < x y z  m n Re Im >  <mR|H|n0>  ==>  < x y z  m n Re Im >  <m0|H|nR>
nk1 = []
nk2 = []
nk3 = []
hop = []
for i in range(0,nat_l):
    for j in range(0,ntot):
        hop.append(float(list[7 + 4 * j + 4 * (ntot+1) * i]))
hop2 = np.zeros(nat_l*ntot)
for i in range(0,9):
    for j in range(nat):
        for k in range(nat):
            for l in range(ntot):
                hop2[l+ntot*(j+nat*(k+nat*i))] = hop[l+ntot*(k+nat*(j+nat*i))]
hop = hop2.copy()
################  change a,b atom from 11 12 13 to 11 21 31... 
################  end change order ################

################   check hermiticity################
hop = set_permutation_symmetry(hop,nat,ntot)
################  end check ermiticity ################
if len(sys.argv) > 1:
    if sys.argv[1] == '1':
        hop = set_simple_asr(hop,nat,ntot)
        print("set asr in simple format") 

################  find wannier_hopping ################
hop3 = np.zeros(nat_l*ntot)
for i in range(0,9):
    for j in range(nat*nat):
            for k in range(ntot):
                hop3[k+ntot*(i+j*9)] = hop2[k+ntot*(j+nat*nat*i)] 
ho = np.array(hop).reshape([nat_l*ntot,1])
sp_poshb = np.zeros([9*nat*nat*ntot,4])
for h in range(9):
    for i in range(nat*nat*ntot):
        for j in range(3):
            sp_poshb[i+nat*nat*ntot*h][j] = sp_posb[i][j]
        sp_poshb[i+nat*nat*ntot*h][3] = ho[i+nat*nat*ntot*h][0]
super_poshb = np.zeros([27*9*nat*nat*ntot,4])
for a in range(27):
    for h in range(9):
        for i in range(nat):
            for j in range(nat*ntot):
                for k in range(3):
                    super_poshb[nat*ntot*(nat*(h+9*a)+i)+j][k] = super_pos[a*nat*ntot+j][k] -cry_pos[j//(ntot)][k]
                super_poshb[nat*ntot*(nat*(h+9*a)+i)+j][3] = ho[j+nat*ntot*(i+h*nat)][0]
count = 0
for i in range(nat*nat*ntot):
    count += multi[i%(nat*ntot)][i//(nat*ntot)]         
atm_ct = np.zeros([nat,nat])
for i in range(nat):
    for j in range(nat):
        for k in range(ntot):
            atm_ct[i][j] += multi[k+ntot*j][i]
mu = np.zeros([count*9,1])
mue = np.zeros([count*9,1])
i = 0
for m in range(nat) :
    for k in range(nat):
        for n in range(ntot):
            mn = multi[n+ntot*k][m]
            for l in range(9):
                mu[i+l*count] = multi[n+ntot*k][m]
                for j in range(mn) :
                    mu[i+j+l*count] = multi[n+ntot*k][m]
            i += mn
for i in range(count*9):
    mue[i] = 1/int(mu[i])
newsp_p = np.zeros([9*count,4])
hr_count=27*ntot
newsphr = np.zeros([9*nat*nat*hr_count,7])
k = 0  

rescale = (1.0E-12)**2*ELECTRONVOLT_SI/(1e-10)**2/AMU_SI/(2*pi)**2     

hmfile=open('./hrmap.txt','w')
hrmap=np.zeros([nat*nat*ntot*nat*27,6])
hmfile.write("{:8d}{:8d}\n".format(int(ntot),int(count)))	

for i in range(nat*nat*ntot):
    sp_dis = sum([c*c for c in sp_bc[i]])
    for j in range(27):
        super_dis = sum([c*c for c in super_bc[i+j*nat*nat*ntot]])
        if abs(super_dis - sp_dis) < 0.00001:
            for l in range(9):
                for m in range(3):
                    newsp_p[k+l*count][m] = super_poshb[i+l*nat*nat*ntot+j*nat_l*ntot][m]
                newsp_p[k+l*count][3] = super_poshb[i+l*nat*nat*ntot+j*nat_l*ntot][3]*mue[k+l*count]
            hrmap[k,0:3]=newsp_p[k,0:3]
            k += 1


kk=0
for i in range(nat):
    for j in range(ntot*nat):
        for k in range(multi[j,i]):
            hrmap[kk,3]=i+1
            hrmap[kk,4]=j//ntot+1
            hrmap[kk,5]=multi[j,i]
            kk +=1
for i in range(count):
    hmfile.write("{:5.0f}{:5.0f}{:5.0f}{:5.0f}{:5.0f}{:5.0f}\n".format(hrmap[i,0],hrmap[i,1],hrmap[i,2],hrmap[i,3],hrmap[i,4],hrmap[i,5]))
hmfile.close()
################  end find wannier_hopping ################

################  generate hr ################
for i in range(3):
        for j in range(3):
            co = 0
            for a in range(nat):
                for b in range(nat):
                    for n in range(int(atm_ct[a][b])):    
                        hrn=int(n+co+count*(j+3*i))    
                        newsp_p[hrn][3] = at_obmass[3*a+j][3*b+i]*(newsp_p[hrn][3])  *(rescale)
                    co += n+1   
hr=np.zeros([9*count,7])
ii = 0
for i in range(3):
    for j in range(3):
        co = 0
        for a in range(nat):
            for b in range(nat):
                for n in range(int(atm_ct[a][b])):    
                    hrn=int(n+co+count*(j+3*i))
                    hr[ii][0:3]=newsp_p[hrn][0:3]
                    hr[ii][3]=3*a+j+1
                    hr[ii][4]=3*b+i+1
                    hr[ii][5]=newsp_p[hrn][3]
                    ii = ii+1
                co += n+1
################  end generate hr ################

################  sorting hopping ################
newsphr[:,0:3] = super_poshb[:,0:3]
newsphr[:,5] = 0
for i in range(27*ntot):
    for j in range(3*nat):
        for k in range(3*nat):
            ii=k+3*nat*(j+3*nat*i)
            newsphr[ii][0]=round(newsphr[ii][0])
            newsphr[ii][1]=round(newsphr[ii][1])
            newsphr[ii][2]=round(newsphr[ii][2])
new_hr1 = np.lexsort([newsphr[:,2],newsphr[:,1],newsphr[:,0]])
newsphr = newsphr[new_hr1,:]
for i in range(27*ntot):
    for j in range(3*nat):
        for k in range(3*nat):
            ii=k+3*nat*(j+3*nat*i)
            newsphr[ii][3]=k+1
            newsphr[ii][4]=j+1
for i in range(27*ntot):
    r=9*nat*nat*i
    for j in range(count):
        if (abs(newsphr[r][0]-hr[j][0])+abs(newsphr[r][1]-hr[j][1])+abs(newsphr[r][2]-hr[j][2])<0.0001):
            for k in range(9):
                for l in range(nat*nat):
                    ii=l+nat*nat*(k+9*i)
                    for kk in range(9):
                        if (abs(newsphr[l+k*nat*nat][3]-hr[j+count*kk][3])+abs(newsphr[l+k*nat*nat][4]-hr[j+count*kk][4])<0.0001):
                            newsphr[ii][5] = hr[j+count*kk][5]
################  end sorting hopping ################                            

################  write wann_hr ################
for scal in range(len(scale[0])):
    file=open('./phhr_{}.dat'.format(scale[1,scal]),'w')
    nrpt=0
    newhr=np.zeros([9*nat*nat*hr_count,7])
    for i in range(27*ntot):
        hr_0 = 0
        for jj in range(9*nat*nat):
            hr_0 = hr_0+abs(newsphr[jj+9*nat*nat*i][5])
        if hr_0 > 1.0E-9:
            for j in range(3*nat):
                for k in range(3*nat):
                    ii=k+3*nat*(j+3*nat*i)
                    newhr[k+3*nat*(j+3*nat*nrpt)] = newsphr[ii]
            nrpt+=1        
    dege_rpts = np.ones(nrpt)
    nl = int(np.ceil(nrpt/15))
    line=" Writen on "+str(datetime.now())+"\n"+"          "+ str(3*nat) + "\n" + "        "+ str(nrpt) + "\n"
    file.write(line)
    for l in range(nl):
        line="    "+'    '.join([str(int(i)) for i in dege_rpts[l*15:(l+1)*15]])
        file.write(line)
        file.write('\n')
    for i in range(nrpt):
        for j in range(3*nat):
            for k in range(3*nat):
                ii=k+3*nat*(j+3*nat*i)
                file.write("{:8.0f}{:8.0f}{:8.0f}{:8.0f}{:8.0f}{:20.10f}{:20.10f}\n".format(newhr[ii][0],newhr[ii][1],newhr[ii][2],newhr[ii][3],newhr[ii][4],scale[0,scal]*newhr[ii][5],scale[0,scal]*newhr[ii][6]))
    file.close()                     
################  end write wann_hr ################


