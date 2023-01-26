CLNSRC=
GINACSRC=

XCFLAGS=${CFLAGS} \
	-O1 -g -Wall -Wextra -std=c++14 -pedantic \
	-pipe -fno-omit-frame-pointer -fpermissive
XLDFLAGS=${LDFLAGS} \
	-lginac -lcln -ldl -lgmp

XCFLAGS_STATIC=${XCFLAGS} -Os -s -static -fdata-sections -ffunction-sections -Wl,--gc-sections
XLDFLAGS_STATIC=${XLDFLAGS}

all: build/fuchsia

build/fuchsia: build/.dir main.cpp fuchsia.cpp Makefile
	date '+static const char VERSION[] = "Fuchsia, built on %Y-%m-%d\n";' >build/version.h
	${CXX} ${XCFLAGS} -include build/version.h -o $@ main.cpp ${XLDFLAGS}

build/fuchsia.static: build/.dir main.cpp fuchsia.cpp Makefile
	env CLN="${CLNSRC}" GINAC="${GINACSRC}" CXX="${CXX}" ./mkversion.sh > build/version_static.h
	${CXX} ${XCFLAGS_STATIC} -include build/version_static.h -o $@ main.cpp ${XLDFLAGS_STATIC}
	upx --best "$@"

test: build/test
	@build/test -a

build/catch.hpp: build/.dir
	wget -O$@ 'https://raw.githubusercontent.com/catchorg/Catch2/v2.5.0/single_include/catch2/catch.hpp'

build/test: build/testrunner.o tests.cpp fuchsia.cpp build/catch.hpp
	${CXX} ${XCFLAGS} -o $@ build/testrunner.o tests.cpp ${XLDFLAGS}

build/testrunner.o: testrunner.cpp build/catch.hpp
	${CXX} ${XCFLAGS} -c -o $@ testrunner.cpp

build/.dir:
	mkdir -p build
	touch $@

clean: phony
	rm -rf build

phony:;
