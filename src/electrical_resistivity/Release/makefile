################################################################################
# Automatically-generated file. Do not edit!
################################################################################


-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include Utilitaries/Visualization/subdir.mk
-include Utilitaries/Storage/subdir.mk
-include Utilitaries/Math/subdir.mk
-include Utilitaries/Geometry/subdir.mk
-include RngStream/subdir.mk
-include InputOutput/subdir.mk
-include FourierModel/subdir.mk
-include EllipseComputation/subdir.mk
-include EPMModel/subdir.mk
-include DomainDefinition/subdir.mk
-include DFNModel/subdir.mk
-include DDPModel/subdir.mk
-include .metadata/.plugins/org.eclipse.cdt.make.core/subdir.mk
-include objects.mk

ifeq ($(MAKECMDGOALS),lib)
-include perform.mk
else ifeq ($(MAKECMDGOALS),Perform)
-include perform.mk
else
-include subdir.mk
endif

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables

.PHONY: exec
exec:
	$(MAKE) BUILD_CXXFLAGS="-O3 -Wall -fopenmp" BUILD_LDFLAGS="-O3 -Wall -fopenmp" ElectricalResistivity

.PHONY: compile
compile:
	$(MAKE) BUILD_CXXFLAGS="-O3 -Wall -fopenmp" BUILD_LDFLAGS="-O3 -Wall -fopenmp" ElectricalResistivity

.PHONY: run
run: ElectricalResistivity
	cd ../.. && ./electrical_resistivity/Release/run_electrical_resistivity.sh

.PHONY: lib
lib:
	$(MAKE) BUILD_CXXFLAGS="-O3 -Wall -fPIC -fopenmp" BUILD_LDFLAGS="-O3 -Wall -fopenmp" Perform

.PHONY: exec-debug
exec-debug:
	$(MAKE) BUILD_CXXFLAGS="-g -Wall -fopenmp" BUILD_LDFLAGS="-g -Wall -fopenmp" ElectricalResistivity

.PHONY: lib-debug
lib-debug:
	$(MAKE) BUILD_CXXFLAGS="-g -Wall -fPIC -fopenmp" BUILD_LDFLAGS="-g -Wall -fopenmp" Perform

.PHONY: test_small_matrix
test_small_matrix:		
	$(MAKE) BUILD_CXXFLAGS="-DTEST_SMALL_MATRIX -g -Wall -fopenmp" BUILD_LDFLAGS="-g -Wall -fopenmp" ElectricalResistivity
	
.PHONY: test_big_matrix
test_big_matrix:		
	$(MAKE) BUILD_CXXFLAGS="-DTEST_BIG_MATRIX -g -Wall -fopenmp" BUILD_LDFLAGS="-g -Wall -fopenmp" ElectricalResistivity

.PHONY: lib_test_small_matrix
lib_test_small_matrix:
	$(MAKE) BUILD_CXXFLAGS="-DLIB_TEST_SMALL_MATRIX -g -Wall -fPIC -fopenmp" BUILD_LDFLAGS="-O3 -Wall -fopenmp" Perform

.PHONY: lib_test_big_matrix
lib_test_big_matrix:
	$(MAKE) BUILD_CXXFLAGS="-DLIB_TEST_BIG_MATRIX -g -Wall -fPIC -fopenmp" BUILD_LDFLAGS="-O3 -Wall -fopenmp" Perform

CXX := g++
PAGURE_ROOT := /home/bcellier/Documents/Stage/new_install/pagure
PAGURE_METIS := $(PAGURE_ROOT)/metis/gcc133/5.1.0
PAGURE_SUITESPARSE := $(PAGURE_ROOT)/suitesparse/gcc133/5.10.1
PAGURE_CGAL := $(PAGURE_ROOT)/cgal/gcc133/5.5.2
PAGURE_BOOST := $(PAGURE_ROOT)/boost/py37/gcc133/1.84.0
PAGURE_PYMOD := $(PAGURE_ROOT)/python-modules/gcc133
PAGURE_PY37 := $(PAGURE_ROOT)/python/gcc133/3.7

PAGURE_INCLUDE_FLAGS := \
	-I$(PAGURE_ROOT)/include \
	-I/home/bcellier/pagure/boost_eigen \
	-I/usr/include/eigen3 \
	-I$(PAGURE_METIS)/include \
	-I$(PAGURE_SUITESPARSE)/include \
	-I$(PAGURE_CGAL)/include \
	-I$(PAGURE_BOOST)/include \
	-I$(PAGURE_PYMOD)/include/python3.7 \
	-I$(PAGURE_PY37)/include/python3.7m

PAGURE_RPATH := $(PAGURE_METIS)/lib:$(PAGURE_SUITESPARSE)/lib:$(PAGURE_CGAL)/lib:$(PAGURE_BOOST)/lib

CXXFLAGS := $(BUILD_CXXFLAGS) $(PAGURE_INCLUDE_FLAGS)
LDFLAGS := $(BUILD_LDFLAGS) -Wl,--no-as-needed -Wl,--copy-dt-needed-entries -Wl,--disable-new-dtags -Wl,-rpath,$(PAGURE_RPATH)

# Tool invocations
ElectricalResistivity: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	$(CXX) $(LDFLAGS) -o ElectricalResistivity $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

Perform: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	$(CXX) $(LDFLAGS) -shared -o "PERFORM.so" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(CC_DEPS)$(C++_DEPS)$(EXECUTABLES)$(OBJS)$(C_UPPER_DEPS)$(CXX_DEPS)$(CPP_DEPS)$(C_DEPS) PERFORM.so ElectricalResistivity
	-@echo ' '

# Test targets
tests : testDFN

testDFN:
	echo "testDFN" > ../../input/SimuType.txt
	rm -f ../../output/testDFN/*
	./ElectricalResistivity
	echo "3frac_SameConduc_test" > ../../input/SimuType.txt
	diff -x .gitkeep -qr ../../TestCases/Output/testDFN/ ../../Output/testDFN/ > /dev/null 2>&1

-include ../makefile.targets
