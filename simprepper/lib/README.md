I copied all this code from here:
https://github.com/JoaoRodrigues/openmm_scripts/tree/master/openmm/lib

I had to make a minor change, though:
I had the change some datatype in the pyx-file from np.float to np.float64

# to make this work

```bash
# this is not included in the environment installation above
micromamba install cython
bash compile_pyx.sh _pwdistance.pyx 
```
