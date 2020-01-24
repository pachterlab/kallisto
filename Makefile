RELEASE_OS ?= local
RELEASE_VERSION ?= local

.PHONY : build install_zlib compile_release_linux compile_release_mac compile_release_windows clean

build:
	# mkdir -p h5ad
	# cd h5ad \
	# && curl -LO https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.15-patch1/src/hdf5-1.8.15-patch1.tar \
	# && tar -xvf hdf5-1.8.15-patch1.tar \
	# && cd hdf5-1.8.15-patch1 \
	# && ./configure --disable-parallel --without-szlib --without-pthread --disable-shared \
	# && make -j \
	# && make install
	- cd ext/htslib \
	&& autoreconf

	mkdir -p build
	cd build \
	&& cmake .. \
	&& make

build_mingw:
	- cd ext/htslib \
	&& autoreconf

	mkdir -p build
	cd build \
	&& cmake .. -DHTSLIB_CONFIGURE_OPTS="--host=x86_64-w64-mingw32" \
	&& make

install_zlib:
	cd ext/zlib \
	&& ./configure --prefix=/usr/src/mxe/usr/x86_64-w64-mingw32.static --static \
	&& make -j \
	&& sudo make install

compile_release_linux compile_release_mac:
	mkdir -p release/kallisto
	cp -rf build/src/kallisto release/kallisto/
	cp -rf license.txt release/kallisto/
	cp -rf README.md release/kallisto/
	cp -rf test release/kallisto/
	cd release \
	&& tar -czvf kallisto_${RELEASE_OS}-${RELEASE_VERSION}.tar.gz kallisto

compile_release_windows:
	mkdir -p release/kallisto
	cp -rf build/src/kallisto.exe release/kallisto/
	cp -rf license.txt release/kallisto/
	cp -rf README.md release/kallisto/
	cp -rf test release/kallisto
	cd release \
	&& zip -r kallisto_${RELEASE_OS}-${RELEASE_VERSION}.zip kallisto

clean:
	rm -rf build
	rm -rf release
	cd ext/htslib && make clean
	cd ext/zlib && make clean
