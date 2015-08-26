---
layout: page
title: "Local build without root"
description: ""
---
{% include JB/setup %}

This tutorial is for building __kallisto__ locally without root access. If you
have root access and a package manager, please see [the other
tutorial](source.html).

#### Requirements:

We will assume you have the following installed:

- A C++11 compiler >= g++-4.8 (might work on g++-4.7, though untested)
    - Since this is very system specific, it is usually best to ask your
      administrator to install this. If not, you should search how to install
      it locally on your system.
- __zlib__ which is installed on most machines.
- __make__ which is also installed on most machines.

### Preliminaries

Since we will be installing things to your home directory, you should add
`$HOME/bin` to your PATH if you haven't already. This ensures that your shell
knows where to look for binaries (cmake, kallisto, etc.). For this current
session, run the following from your terminal:

```
export PATH=$HOME/bin:$PATH
```

Afterwards, place the same code into your shell startup file (e.g. one of ~/.bashrc, ~/.zshrc, etc.).

All downloads can go somewhere, let's say `~/downloads`.

### Building and installing CMake

The easiest way to install CMake is from source. Head over to the [CMake
downloads page](http://www.cmake.org/download/) and get the latest "Unix/Linux
Source" \*.tar.gz file. Make sure to set the `--prefix` flag correctly,
otherwise you won't have permissions to install files at the default location.

`tar -xf cmake*.tar.gz`

`cd cmake*`

`./configure --prefix=$HOME`

`make`

`make install`

You should now have the most up-to-date installation of cmake. Check the
version by typing:

```
cmake --version
```

### Building and installing HDF5

Download 'configure' version of [HDF5](https://www.hdfgroup.org/HDF5/release/obtainsrc.html#conf):

`wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.15-patch1.tar.gz`

Extract the tarball and move to the directory:

`tar -xf hdf5-1.8.*`

`cd hdf5-1.8.*`

Now, run `configure` with the following options:

`./configure --disable-parallel --without-szlib --without-pthread --prefix=$HOME`

Compile and install:

`make`

`make install`

All of the important HDF5 tools will be at `$HOME/bin` and libraries/include
files at: `$HOME/lib` and `$HOME/include`.

### Building and installing kallisto

Now that we have the depencies installed, __kallisto__ should install pretty
easily. Download the latest source tarball [from
GitHub](https://github.com/pachterlab/kallisto/releases), extract, and change
directories:

`wget ...` (this is the GitHub 'Source code (tar.gz)' link)

`tar -xf v0.4*`

`cd kallisto-*`

Make a build directory and move to it

`mkdir build`

`cd build`

Since we put everything in `$HOME`, CMake is smart enough to look there and the
following should work pretty seamlessly. Make sure you run CMake with the
following command to install the binaries in `$HOME`:

`cmake -DCMAKE_INSTALL_PREFIX=$HOME ..`

Next, build and install:

`make`

`make install`

You will now have `kallisto` in `$HOME/bin/kallisto`.
