# IR2PW and IR2TR
to compute Irreducible Representations by linking the IRVSP library (libIRVSP.a) to PW code and TB model.</br>


# libIRVSP.a 
The library is created by IRVSP (v2), according to the CRTs on the Bilbao Crystallographic Server (BCS). </br>
The library can be linked to by DFT packages, such as Quantum Espresso, VASP, Siesta, Abinit, ELK and Wien2k etc.

* lib_irrep_bcs.tar.gz : The IRVSP library released on 2023.x.XX.
* how to make:

      $$   sdfasdf
      $$   asfdasdf
      $$  sdfasdf

# IR2PW (with interface to QE, VASP)
to compute irreducible representations with interface to plane-wave (PW) DFT packages.

* src_ir2pw_QE.tar.gz  :interface to the Quantum-Espresso package.</br>
required files: nscf_b.out and $outdir (output directory)
                      

* src_ir2pw_VASP.tar.gz:interface to the Vienna ab-initio Simulation Package.</br>
required files: OUTCAR and WAVECAR


# IR2TB (with interface to Wannier90，QE/ph.x)
to compute irreducible representations with interface to orthogonal tight-binding (TB) models (Wannier90, Slater-Koster, Phonon TB). </br>
It works for phonon, elctron, magnon systems.

* src_ir2tb_wann.tar.gz : interface to orthogonal TB model. </br>
required files: tbbox.in and wann_hr.dat </br>
axiliaries: 

          pwscf2tbbox.py : to convert pwscf.out (QE) to tbbox.in 
          fc2hr.py       : to convert  xx.fx to wann_hr.dat

* src_ir2tb_phx.tar.gz : interface to QE phonon (q) calculations. </br>
required files: ph_q.out and $fildyn (ph.dyn0, ph.dyn1, ph.dyn...) </br>
axiliaries: 

          dyn2wf.py  : to convert  .dyn* to phonon wavefunctions
          mode2wf.py : to convert .modes to phonon wavefunctions 
