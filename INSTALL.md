# kallisto installation

Requirements
------------
- CMake version >= 2.8.12
    - Can be installed via homebrew: `brew install cmake`
- zlib (should be installed on OSX >= 10.9)
- HDF5 C libraries

Installation
------------

1. First clone the repository:

    ```
    git clone https://github.com/pachterlab/kallisto.git
    ```

1. Move to the source directory:

    ```
    cd kallisto
    ```

1. Make a build directory and move there

    ```
    mkdir build
    cd build
    ```

1. Execute cmake. There are a few options:
    - `-DCMAKE_INSTALL_PREFIX:PATH=$HOME` which will put kallisto in
       `$HOME/bin` as opposed to the default (`/usr/local/bin`)

    ```
    cmake ..
    ```

    This is only required when one of the `CMakeLists.txt` files changes or new
    source files are introduced. It will make a new set of `Makefile`s

1. Build the code

    ```
    make
    ```

    Optionally install into the cmake install prefix path:

    ```
    make install
    ```
1. Run the code. The source tree hierarchy is copied, so you can find kallisto
   at `src/kallisto` and tests at `test/tests`.

After performing these steps, you can simply build using make as long as new
source files aren't introduced or the `CMakeLists.txt` scripts aren't modified.
