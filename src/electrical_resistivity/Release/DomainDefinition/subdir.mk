################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../DomainDefinition/Domain.cpp \
../DomainDefinition/ElectricProperties.cpp 

OBJS += \
./DomainDefinition/Domain.o \
./DomainDefinition/ElectricProperties.o 

CPP_DEPS += \
./DomainDefinition/Domain.d \
./DomainDefinition/ElectricProperties.d 


# Each subdirectory must supply rules for building sources it contributes
DomainDefinition/%.o: ../DomainDefinition/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(CXXFLAGS) -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


