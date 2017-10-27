XCFLAGS=${CFLAGS} -O0 -g -Wall -Wextra -std=c++14 -pedantic -fno-omit-frame-pointer
XLDFLAGS=${LDFLAGS} -lginac -lcln -ldl

all: _build/main

run: _build/main
	@_build/main

_build/main: _build/.dir main.cpp fuchsia.cpp Makefile
	${CXX} ${XCFLAGS} -o $@ main.cpp ${XLDFLAGS}

test: _build/test
	@_build/test

_build/test: _build/testrunner.o tests.cpp fuchsia.cpp
	${CXX} ${XCFLAGS} -o $@ _build/testrunner.o tests.cpp ${XLDFLAGS}

_build/testrunner.o: _build/.dir testrunner.cpp
	${CXX} ${XCFLAGS} -c -o $@ testrunner.cpp

_build/.dir:
	ln -sf "$$(mktemp -d)" _build
	touch $@

clean: phony
	rm -rf $$(readlink -f _build) _build

phony:;
