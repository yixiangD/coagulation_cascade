# Coagulation Cascade
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
>> ../../lmp_stam < chan.equil.run
```
### 2nd: add a driven force and allow chemical species propagattion
```
>> ../../lmp_stam < chan.rheo.run
```
### 3nd: turn on adhesion force, enable platelet adhesion on wall
```
>> ../../lmp_stam < chan.adr.run
```

# Copyright

**You are not allowed to use this repo for any commercial use or academic
publication, unless a permission from the repo owner is given.**


**If you use this repo for your academic research, we require a citation to the following publication.**

TBA
