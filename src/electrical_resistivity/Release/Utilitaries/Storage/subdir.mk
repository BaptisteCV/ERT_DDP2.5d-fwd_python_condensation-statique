################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Utilitaries/Storage/Structures.cpp \
../Utilitaries/Storage/UblasStructures.cpp \
../Utilitaries/Storage/UtilitaryStorage.cpp 

OBJS += \
./Utilitaries/Storage/Structures.o \
./Utilitaries/Storage/UblasStructures.o \
./Utilitaries/Storage/UtilitaryStorage.o 

CPP_DEPS += \
./Utilitaries/Storage/Structures.d \
./Utilitaries/Storage/UblasStructures.d \
./Utilitaries/Storage/UtilitaryStorage.d 


# Each subdirectory must supply rules for building sources it contributes
Utilitaries/Storage/%.o: ../Utilitaries/Storage/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(CXXFLAGS) -I /home/bcellier/Documents/Stage/new_install/pagure/include -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


