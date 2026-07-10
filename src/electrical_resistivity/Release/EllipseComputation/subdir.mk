################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../EllipseComputation/EllipseComputation.cpp 

OBJS += \
./EllipseComputation/EllipseComputation.o 

CPP_DEPS += \
./EllipseComputation/EllipseComputation.d 


# Each subdirectory must supply rules for building sources it contributes
EllipseComputation/%.o: ../EllipseComputation/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(CXXFLAGS) -I /home/bcellier/Documents/Stage/new_install/pagure/include -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


