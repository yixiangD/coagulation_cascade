# LAMMPS for Coagulation Cascade Modeling
## Compile source code
```
>> cd src/
>> make
```
## Running simulations
```
>> cd ../case/nrbc # normal case, change to drbc for diabetic case
```
### 1st: run the system to equilibrium, without force
```
>> ../../src/lmp_stam < chan.equil.run
```
### 2nd: add a driven force and allow chemical species advection and diffusion
```
>> ../../src/lmp_stam < chan.rheo.run
```
### 3nd: turn on adhesion force, enable platelet adhesion on wall
```
>> ../../src/lmp_stam < chan.adr.run
```

# Copyright

**You are not allowed to use this repo for any commercial use or academic
publication, unless a permission from the repo owner is given.**


**If you use this repo for your academic research, we require a citation to the following publication.**
```
@article{yazdani2021integrating,
  title={Integrating blood cell mechanics, platelet adhesive dynamics and coagulation cascade for modelling thrombus formation in normal and diabetic blood},
  author={Yazdani, Alireza and Deng, Yixiang and Li, He and Javadi, Elahe and Li, Zhen and Jamali, Safa and Lin, Chensen and Humphrey, Jay D and Mantzoros, Christos S and Em Karniadakis, George},
  journal={Journal of the Royal Society Interface},
  volume={18},
  number={175},
  pages={20200834},
  year={2021},
  publisher={The Royal Society}
```
