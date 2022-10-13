# IR2PW and IR2TB
to compute irreducible representations by linking the library (libIRVSP.a) to PW code and TB model.</br>


# libIRVSP.a 
* The library is created by IRVSP, according to the CRTs on the Bilbao Crystallographic Server (BCS).
* The library can be linked to by DFT packages, such as Quantum Espresso, VASP, Siesta, Abinit, ELK and Wien2k etc.

# ir2pw (with interface to QE, VASP)
to compute irreducible representations with interface to plane-wave (PW) DFT packages.

src_ir2pw_QE.tar.gz  :interface to the Quantum-Espresso package.</br>
required files: nscf_b.out and "output" directory
                      

src_ir2pw_VASP.tar.gz:interface to the Vienna ab-initio Simulation Package.</br>
required files: OUTCAR and WAVECAR


# ir2tb
to compute irreducible representations with interface to tight-binding (TB) models. </br>
It works for phonon, elctron, magnon systems.

src_ir2tb.tar.gz : interface to Wannier90-based TB model.
                     
