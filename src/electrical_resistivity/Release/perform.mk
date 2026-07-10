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
	$(CXX) $(CXXFLAGS) -fPIC -I /home/bcellier/Documents/Stage/new_install/pagure/include -I /home/bcellier/Documents/Stage/new_install/pagure/python-modules/gcc133/include/python3.7 -I /home/bcellier/Documents/Stage/new_install/pagure/python-modules/gcc133/include/python3.7/pybind11 -I /home/bcellier/Documents/Stage/new_install/pagure/python/gcc133/3.7/include/python3.7m -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '
