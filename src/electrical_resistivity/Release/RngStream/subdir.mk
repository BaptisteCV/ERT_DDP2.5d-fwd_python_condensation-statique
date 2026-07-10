################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../RngStream/rngstream.cpp 

OBJS += \
./RngStream/rngstream.o 

CPP_DEPS += \
./RngStream/rngstream.d 


# Each subdirectory must supply rules for building sources it contributes
RngStream/%.o: ../RngStream/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(CXXFLAGS)  -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


