# kallisto

Requirements
------------
- CMake version >= 2.8.8
    - Can be installed via homebrew: `brew install cmake`
- zlib (should be installed on OSX >= 10.9)

Installation
------------

1. First clone the repository:

    ```
    git clone git@github.com:pimentel/kallisto.git
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

1. Execute cmake

    ```
    cmake ..
    ```

    This is only required when one of the `CMakeLists.txt` files changes or new
    source files are introduced. It will make a new set of `Makefile`s
1. Build the code

    ```
    make
    ```
1. Run the code. The source tree hierarchy is copied, so you can find kallisto
   at `src/kallisto` and tests at `test/tests`.

After performing these steps, you can simply build using make as long as new
source files aren't introduced or the `CMakeLists.txt` scripts aren't modified.
