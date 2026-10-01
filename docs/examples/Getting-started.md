---
title: Getting started
---


## Compilation of the C++ code


### Basics

The C++ code is built using rules given by cmake. In order to compile the code on your computer, you'll therefore need to have cmake installed, as well as a C++ compiler. 

### Dependencies

The main dependency is Eigen. Version 5.0.0 of Eigen is used, as it is the newest version in 2026, when the majority of the code was written. Some breaking changes were introduced in version 5.0.0 such that if you have version 3.4 installed on your computer, the code won't compile. Also, note that the version naming changed in the resease (see [notes](https://gitlab.com/libeigen/eigen/-/releases/5.0.0)), such that it is normal that Eigen is still imported as Eigen3 in the cmake (ignore this last phrase if you did'nt bother about that before). You can install it from your package manager or git clone the repo and include the include dir (Eigen is a header only library, so no need to build it).

If you wish to use result storage, you also need to have hdf5 installed on your system. Depending on your OS:

-   Macos: `brew install hdf5` should do the trick.
-   Linux: something like `apt install hdf5-tools` should work if you have apt as a package manager.
-   Windows: no idea.

I fyou don't want to use result storage, then you can still build the rest of the code without having to install hdf5 by setting the TARTE_USE_HDF5 flag to OFF in the root CMakeLists.txt:

`option(TARTE_USE_HDF5 "Enable HDF5 support" OFF)`

### Command line compilation

Running the following command lines in a bash-like shell should build the code:

```
mkdir build 
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

The first two lines create a build folder and move into it, then cmake generates the building rules (here in release mode), and make builds the executables.


### Running examples executable 

After building, examples can be run as executables. 

## Python configuration 

Some python code is provided to simplify generation and analysis of simulation data. All of the python code relies on the C++ part being built with hdf5 support on, such that the previous steps are needed before trying to run the python part. Once this is out of the way, and assuming you have python and a package manager such as pip installed on your computer, running the following commands form the root of the repository will create a virtual environment with the necessary dependencies installed:
```
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

You can check that everything is working by then trying to run a simulation from an example python script:
```
cd python
python VoiceReadWrite.py
```
which should generate a `ReadWriteExample.hdf5` file in the python directory with voice simulation results in it.

## Building the docs

If for some reason you wish to build and serve the present documentation from your computer, after activating the python virtual environment as described in the previous section, you can run 
```
mkdocs serve
```
from the root of the repository.