RELEASE_OS ?= local
RELEASE_VERSION ?= local

.PHONY : build build_mingw install_zlib compile_release_linux compile_release_mac compile_release_windows clean

build:
	mkdir -p build
	cd build \
	&& cmake .. \
	&& make -j

build_mingw:
	mkdir -p build
	cd build \
	&& cmake .. -DHTSLIB_CONFIGURE_OPTS="--host=x86_64-w64-mingw32" \
	&& make -j

install_zlib:
	mkdir -p ext/zlib \
	&& cd ext/zlib \
	&& wget https://zlib.net/fossils/zlib-1.2.5.tar.gz \
	&& tar -xvf zlib-1.2.5.tar.gz --strip 1
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
