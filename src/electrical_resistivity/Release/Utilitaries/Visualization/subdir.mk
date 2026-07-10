################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Utilitaries/Visualization/DisplayResults.cpp \
../Utilitaries/Visualization/VisuStructures.cpp \
../Utilitaries/Visualization/VisuTools.cpp 

OBJS += \
./Utilitaries/Visualization/DisplayResults.o \
./Utilitaries/Visualization/VisuStructures.o \
./Utilitaries/Visualization/VisuTools.o 

CPP_DEPS += \
./Utilitaries/Visualization/DisplayResults.d \
./Utilitaries/Visualization/VisuStructures.d \
./Utilitaries/Visualization/VisuTools.d 


# Each subdirectory must supply rules for building sources it contributes
Utilitaries/Visualization/%.o: ../Utilitaries/Visualization/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(CXXFLAGS) -I /usr/include/GL/ -I /home/bcellier/Documents/Stage/new_install/pagure/include -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


