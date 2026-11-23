################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables
CPP_SRCS += \
../PERFORM.cpp

OBJS += \
./PERFORM.o

CPP_DEPS += \
./PERFORM.d


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(CXXFLAGS) -fPIC $(EXTRA_INCLUDE_FLAGS) -I $(PYMOD_INCLUDE) -I $(PYMOD_INCLUDE)/pybind11 -I $(PYTHON_INCLUDE) -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '
