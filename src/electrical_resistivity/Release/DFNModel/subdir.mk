################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../DFNModel/DFNComputation.cpp \
../DFNModel/DFNGeneration.cpp \
../DFNModel/FractureMesh.cpp \
../DFNModel/NetworkMeshes.cpp 

OBJS += \
./DFNModel/DFNComputation.o \
./DFNModel/DFNGeneration.o \
./DFNModel/FractureMesh.o \
./DFNModel/NetworkMeshes.o 

CPP_DEPS += \
./DFNModel/DFNComputation.d \
./DFNModel/DFNGeneration.d \
./DFNModel/FractureMesh.d \
./DFNModel/NetworkMeshes.d 


# Each subdirectory must supply rules for building sources it contributes
DFNModel/%.o: ../DFNModel/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(CXXFLAGS) -I /home/bcellier/Documents/Stage/new_install/pagure/include -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


