# IR2PW and IR2TB
to compute Irreducible Representations by linking the IRVSP library (libIRVSP.a) to PW code and TB model.</br>


# libIRVSP.a 
This library is created by IRVSP (v2), according to the CRTs on the Bilbao Crystallographic Server (BCS). </br>
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

* "ir2pw" in src_ir2pw_QE.tar.gz: interface to the Quantum-Espresso package.</br>
required files: nscf_b.out and $outdir (output directory) </br>
 \$ set nscf_b.in; ph.x -in nscf_b.in >nscf_b.out

      calculation = 'bands'
       outdir     = './tmp'
      verbosity   = 'high'


* "irvsp" in src_ir2pw_VASP.tar.gz: interface to the Vienna ab-initio Simulation Package.</br>
required files: OUTCAR and WAVECAR


# IR2TB (with interface to Wannier90，QE/ph.x)
to compute irreducible representations with interface to orthogonal tight-binding (TB) models (Wannier90, Slater-Koster, Phonon TB). </br>
It works for phonon, electron, magnon systems.

* "ir2tb" in src_ir2tb_hr.tar.gz : interface to orthogonal TB model. </br>
required files: tbbox.in and lda_hr.dat/soc_hr.dat </br>
axiliaries: pwscf2tbbox.py and fc2hr.py

          pwscf2tbbox.py : to convert pwscf.out (QE) to tbbox.in 
          fc2hr.py       : to convert ph.fc to phonon TB ph_hr.dat

* "ir2ph" in src_ir2tb_ph.tar.gz : interface to QE phonon (q) calculations. </br>
required files: tbbox.in and ph_wf.dat </br>
axiliaries: IRphx.sh <br>
 \$ vim ph.inp <br>
 \$ IRph.sh <br>
 \$ ph.x < q*.inp > q*-out <br>
 \$ IRphx.sh 
 
          IRphx.sh : ph.inp and/or q*.dyn/q*-out
          1. to prepare q*.inp (from ph.inp)
          2. to convert q*.dyn to ph_wf.dat


# Notices
* Please put all these folders in the same directory in order to make it successfully.
* Make the IRVSP library first. Then make IR2PW/IR2TB.
* ir2pw/irvsp/ir2tb/ir2ph -sg xx -nb xx xx >outir2
