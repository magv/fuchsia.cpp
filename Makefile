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

build/fuchsia: build/.dir main.cpp fuchsia.cpp Makefile version.h
	${CXX} ${XCFLAGS} -DVERSION_H=\"version.h\" -o $@ main.cpp ${XLDFLAGS}

build/fuchsia.static: build/.dir main.cpp fuchsia.cpp Makefile version_static.h
	${CXX} ${XCFLAGS_STATIC} -DVERSION_H=\"version_static.h\" -o $@ main.cpp ${XLDFLAGS_STATIC}
	upx --best "$@"

test: build/test
	@build/test

build/test: build/testrunner.o tests.cpp fuchsia.cpp
	${CXX} ${XCFLAGS} -o $@ build/testrunner.o tests.cpp ${XLDFLAGS}

build/testrunner.o: build/.dir testrunner.cpp
	${CXX} ${XCFLAGS} -c -o $@ testrunner.cpp

build/.dir:
	ln -sf "$$(mktemp -d)" build
	touch $@

version.h: phony
	date '+"Fuchsia, built on %Y-%m-%d\n"' > $@

version_static.h: phony
	printf "%s" 'R"(' > $@
	env CLN="${CLNSRC}" GINAC="${GINACSRC}" ./mkversion.sh >> $@
	printf "%s" ')"' >> $@

clean: phony
	rm -rf $$(readlink -f build) build version.h version_static.h

phony:;
