CLNSRC=
GINACSRC=

XCFLAGS=${CFLAGS} \
	-O0 -g -Wall -Wextra -std=c++14 -pedantic \
	-pipe -fno-omit-frame-pointer -fpermissive \
	-fdata-sections -ffunction-sections -Wl,--gc-sections
XLDFLAGS=${LDFLAGS} \
	-lginac -lcln -ldl -lgmp

XCFLAGS_STATIC=${XCFLAGS} -Os -s -static
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

build/test: build/testrunner.o tests.cpp fuchsia.cpp
	${CXX} ${XCFLAGS} -o $@ build/testrunner.o tests.cpp ${XLDFLAGS}

build/testrunner.o: build/.dir testrunner.cpp
	${CXX} ${XCFLAGS} -c -o $@ testrunner.cpp

build/.dir:
	ln -sf "$$(mktemp -d)" build
	touch $@

clean: phony
	rm -rf $$(readlink -f build) build

phony:;
