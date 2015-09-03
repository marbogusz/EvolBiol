# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/indestimate.cpp 

OBJS += \
./src/indestimate.o 

CPP_DEPS += \
./src/indestimate.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	$(CC) $(INC_PAR) $(CPPFLAGS) -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


