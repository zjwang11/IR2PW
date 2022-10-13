# IR2PW and IR2TB
to compute irreducible representations by linking the library (libIRVSP.a) to PW code and TB model.</br>


# libIRVSP.a 
The library is created by IRVSP, according to the CRTs on the Bilbao Crystallographic Server (BCS). </br>
The library can be linked to by DFT packages, such as Quantum Espresso, VASP, Siesta, Abinit, ELK and Wien2k etc.

* lib_irrep_bcs.tar.gz : The IRVSP library released on 2023.x.XX.

# ir2pw (with interface to QE, VASP)
to compute irreducible representations with interface to plane-wave (PW) DFT packages.

* src_ir2pw_QE.tar.gz  :interface to the Quantum-Espresso package.</br>
required files: nscf_b.out and "output" directory
                      

* src_ir2pw_VASP.tar.gz:interface to the Vienna ab-initio Simulation Package.</br>
required files: OUTCAR and WAVECAR


# ir2tb
to compute irreducible representations with interface to orthogonal tight-binding (TB) models. </br>
It works for phonon, elctron, magnon systems.

* src_ir2tb_ort.tar.gz : interface to orthogonal TB model.</br>
required files: tbbox.in and hr.dat
                     
* src_ir2tb_ort.tar.gz : interface to QE phonon q calculations. </br>
required files: ph_q.out and dynamic_x.xxx
