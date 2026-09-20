SHELL := /bin/sh

# Compiler discovery for galaxy-parallel native workers.
DEPENDENCY_ROOT ?= .
GCC_PACKAGE := $(firstword $(wildcard $(DEPENDENCY_ROOT)/.conda-pkgs/gcc_impl_osx-arm64-*/bin/arm64-apple-darwin*-gcc))
GCC_LIBDIR := $(firstword $(wildcard $(DEPENDENCY_ROOT)/.conda-pkgs/libgcc-devel_osx-arm64-*/lib/gcc/arm64-apple-darwin*/*))
CC := $(if $(GCC_PACKAGE),$(GCC_PACKAGE),gcc)

CPPFLAGS := -I$(DEPENDENCY_ROOT)/.deps/include -I$(DEPENDENCY_ROOT)/vendor/cubature
CFLAGS := -std=c11 -D_DEFAULT_SOURCE -O3 -fopenmp -Wall -Wextra \
          -Wno-unused-function -Wno-unused-parameter
LDLIBS := $(DEPENDENCY_ROOT)/.deps/lib/libgsl.a $(DEPENDENCY_ROOT)/.deps/lib/libgslcblas.a $(DEPENDENCY_ROOT)/.deps/lib/libcubature.a \
          -L$(GCC_LIBDIR) -L$(DEPENDENCY_ROOT)/.deps/lib \
          -lemutls_w -lm -Wl,-rpath,$(abspath $(DEPENDENCY_ROOT))/.deps/lib

.DEFAULT_GOAL := all
.PHONY: all
all: solver

# Build the interface only; reference executables are maintained externally.
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
LIB_SUFFIX := dylib
SHARED_FLAGS := -dynamiclib
else
LIB_SUFFIX := so
SHARED_FLAGS := -shared
endif
PYTHON ?= python3
.PHONY: shared
# Compatibility build target; there is no separate preparation library.
shared: solver

build:
	mkdir -p $@

.PHONY: test-portable test
test-portable: solver
	PYTHONPATH=src $(PYTHON) -m unittest discover -s tests -p test_preparation.py -v

test: solver test-runtime
	PYTHONPATH=src $(PYTHON) -m unittest discover -s tests -v

.PHONY: solver
# Build only the worker library; no reference executable is required.
# All native physics dependencies live under csrc/physics.
SOLVER_LIBRARY := build/libcongruents_solver.$(LIB_SUFFIX)
ifeq ($(UNAME_S),Darwin)
# Keep GSL implementation symbols private to avoid sharing its global error
# handler with unrelated libraries loaded into the same Python process.
SOLVER_LINK_FLAGS := -Wl,-dead_strip -Wl,-exported_symbol,_cg_solver_abi -Wl,-exported_symbol,_cg_solver_openmp_enabled -Wl,-exported_symbol,_cg_properties -Wl,-exported_symbol,_cg_transport -Wl,-exported_symbol,_cg_spectra
else
SOLVER_LINK_FLAGS := -Wl,--gc-sections -Wl,--exclude-libs,ALL
endif
solver: $(SOLVER_LIBRARY)

.PHONY: test-runtime test-runtime-sanitized
build/test_solver_runtime: tests/test_solver_runtime.c csrc/solver.c csrc/solver.h csrc/solver_runtime.h $(wildcard csrc/physics/*.h) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LDLIBS) -o $@

test-runtime: build/test_solver_runtime
	build/test_solver_runtime

# Audit the same ABI serially without requiring a sanitizer-compatible OpenMP runtime.
SANITIZER_CC ?= clang
ifeq ($(UNAME_S),Darwin)
RUNTIME_ASAN_OPTIONS ?= detect_leaks=0
else
RUNTIME_ASAN_OPTIONS ?= detect_leaks=1
endif
build/test_solver_runtime_sanitized: tests/test_solver_runtime.c csrc/solver.c csrc/solver.h csrc/solver_runtime.h $(wildcard csrc/physics/*.h) | build
	$(SANITIZER_CC) $(CPPFLAGS) -std=c11 -D_DEFAULT_SOURCE -DCG_SERIAL_AUDIT -Wno-unknown-pragmas -g -O1 -fno-omit-frame-pointer -fsanitize=address,undefined $< $(LDLIBS) -o $@

test-runtime-sanitized: build/test_solver_runtime_sanitized
	ASAN_OPTIONS=$(RUNTIME_ASAN_OPTIONS) build/test_solver_runtime_sanitized

.PHONY: test-runtime-leaks
test-runtime-leaks: build/test_solver_runtime
	leaks --atExit -- build/test_solver_runtime

$(SOLVER_LIBRARY): Makefile csrc/solver.c csrc/solver.h csrc/solver_runtime.h $(wildcard csrc/physics/*.h) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) -ffunction-sections -fdata-sections -fPIC -fvisibility=hidden $(SHARED_FLAGS) csrc/solver.c $(LDLIBS) $(SOLVER_LINK_FLAGS) -o $@

# Optional numerical comparison; neither normal builds nor runs invoke it.
REFERENCE_OUTPUT ?= ../CONGRUENTS1.5/output
PYTHON_OUTPUT ?= output/production
.PHONY: compare
compare:
	$(PYTHON) tests/compare_outputs.py "$(REFERENCE_OUTPUT)" "$(PYTHON_OUTPUT)"
