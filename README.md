# IR2PW and IR2TB
to compute Irreducible Representations by linking the IRVSP library (libIRVSP.a) to PW code and TB model.</br>


# libIRVSP.a 
The library is created by IRVSP (v2), according to the CRTs on the Bilbao Crystallographic Server (BCS). </br>
It can be linked to by DFT packages, such as Quantum Espresso, VASP, Siesta, Abinit, ELK and Wien2k etc. </br>
See the library details in Ref: J. Gao, et al. Comput. Phys. Comm. 261, 107760 (2021). https://doi.org/10.1016/j.cpc.2020.107760.

* lib_irrep_bcs.tar.gz : The IRVSP library is released on 2022.10.18.
* how to make:

      $  tar -zxvf lib_irrep_bcs.tar.gz
      $  cd lib_irrep_bcs
      $  ./configure.sh
      $  make

# IR2PW (with interface to QE, VASP)
to compute irreducible representations with interface to plane-wave (PW) DFT packages.

* src_ir2pw_QE.tar.gz: interface to the Quantum-Espresso package.</br>
required files: nscf_b.out and $outdir (output directory) </br>
 \$ set nscf_b.in; ph.x -in nscf_b.in >nscf_b.out

      calculation = 'bands'
       outdir     = './tmp'
      verbosity   = 'high'


* src_ir2pw_VASP.tar.gz: interface to the Vienna ab-initio Simulation Package.</br>
required files: OUTCAR and WAVECAR


# IR2TB (with interface to Wannier90，QE/ph.x)
to compute irreducible representations with interface to orthogonal tight-binding (TB) models (Wannier90, Slater-Koster, Phonon TB). </br>
It works for phonon, electron, magnon systems.

* src_ir2tb_hr.tar.gz : interface to orthogonal TB model. </br>
required files: tbbox.in and ldawann_hr.dat </br>
axiliaries: 

          pwscf2tbbox.py : to convert pwscf.out (QE) to tbbox.in 
          fc2hr.py       : to convert ph.fc to phonon TB wann_hr.dat

* src_ir2tb_ph.tar.gz : interface to QE phonon (q) calculations. </br>
required files: q*.out and q*.dyn (q0.dyn, q1.dyn, q2.dyn...) </br>
axiliaries: 

          dyn2wf.py  : to convert  q*.dyn to phonon TB wavefunctions


# Notices
* Please put all these folders in the same directory in order to make it successfully.
* Make the IRVSP library first. Then make IR2PW/IR2TB.
* ir2pw/irvsp/ir2tb/ir2ph -sg xx -nb xx xx >outir2
