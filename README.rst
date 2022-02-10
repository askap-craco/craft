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




For Development
##############
```
git clone ...
cd craft
python3.8 -v venv venv
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
