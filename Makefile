XCFLAGS=${CFLAGS} -O0 -g -Wall -Wextra -std=c++14 -pedantic -fno-omit-frame-pointer
XLDFLAGS=${LDFLAGS} -lginac -lcln -ldl

all: build/main

run: build/main
	@build/main

build/main: build/.dir main.cpp fuchsia.cpp Makefile
	${CXX} ${XCFLAGS} -o $@ main.cpp ${XLDFLAGS}

test: build/test
	@build/test

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
