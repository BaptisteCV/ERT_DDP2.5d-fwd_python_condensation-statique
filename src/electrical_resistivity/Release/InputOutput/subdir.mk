################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../InputOutput/BoundaryConditions.cpp \
../InputOutput/Parameters.cpp \
../InputOutput/Results.cpp 

OBJS += \
./InputOutput/BoundaryConditions.o \
./InputOutput/Parameters.o \
./InputOutput/Results.o 

CPP_DEPS += \
./InputOutput/BoundaryConditions.d \
./InputOutput/Parameters.d \
./InputOutput/Results.d 


# Each subdirectory must supply rules for building sources it contributes
InputOutput/%.o: ../InputOutput/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(CXXFLAGS) -I /home/bcellier/Documents/Stage/new_install/pagure/include -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


