CRAFT
##########

This project has utilities for CRACO project

``pip install --upgrade .``


Dependencies
############
This uses numba which requires llvmlite. The LLLVMLITE library requires you have a particular version of llvm installed on your system. See for a mapping https://pypi.org/project/llvmlite/

I have LLVM 6.0 on my machine. Then I can install the rest of the the stuff.

```
pip install llvmlite==0.26.0
```

But this didn't work for me.

In the end I followed these instructions to get LLVM11 installed on my Ubuntu 18.04LTS box https://gist.github.com/kittywhiskers/a3395cb41206d8aa777ce0a8b722d37e then installed numba with

```
LLVM_CONFIG=llvm-config-11 pip install numba

```

On python3.6 this worked for me
###############################
```
sudo apt install llvm-10
LLVM_CONFIG=llvm-config-10 pip install numba
```

LAPACK etc
######
```
sudo apt install libblas3 liblapack3 liblapack-dev libblas-dev gfortran libatlas-base-dev python-numpy python-scipy

```

To run the pipeline
###################
Creating test data currently requires ATNF Mirid. Install it and make sure it's in your path. See https://www.atnf.csiro.au/computing/software/miriad/INSTALL.html


For Development
##############
```
git clone ...
cd craft
# pick up system scipy otherwise it takes ages to install
python3.8 -v venv venv  --system-site-packages
source venv/bin/activate

# the following pip installs all teh required dependencies and installs
# craft in editable mode that means you can edit scripts in the craft
# directory and you don't have to re-instally it to see those changes.
pip install -e . 
```

# For notebooks

You probably want to install jupyter

```
pip install jupyter
cd notebooks
jupyter notebook
```
