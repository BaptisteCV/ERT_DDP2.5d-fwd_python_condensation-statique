################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../EPMModel/EPMComputation.cpp \
../EPMModel/EPMComputationCSEigen.cpp \
../EPMModel/EPMDefinition.cpp \
../EPMModel/EPMDefinitionCSEigen.cpp 

OBJS += \
./EPMModel/EPMComputation.o \
./EPMModel/EPMComputationCSEigen.o \
./EPMModel/EPMDefinition.o \
./EPMModel/EPMDefinitionCSEigen.o 

CPP_DEPS += \
./EPMModel/EPMComputation.d \
./EPMModel/EPMComputationCSEigen.d \
./EPMModel/EPMDefinition.d \
./EPMModel/EPMDefinitionCSEigen.d 


# Each subdirectory must supply rules for building sources it contributes
EPMModel/%.o: ../EPMModel/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(CXXFLAGS) -I /home/bcellier/Documents/Stage/new_install/pagure/include -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


