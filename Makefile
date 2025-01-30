# Basic dvm-dos-tem Makefile 

# Add compiler flag for enabling floating point exceptions:
# -DBSD_FPE for BSD (OSX)
# -DGNU_FPE for various Linux

CC=g++
CFLAGS=-c -ansi -g -gdwarf-2 -std=c++11 -fPIC -DBOOST_ALL_DYN_LINK -Werror # -W -Wall -Werror -Wno-system-headers
LIBS=-lnetcdf -lboost_system -lboost_filesystem \
-lboost_program_options -lboost_thread -lboost_log -ljsoncpp -lpthread -lreadline -llapacke

USEMPI = false
USEOMP = false

ifeq ($(USEMPI),true)
  MPIINCLUDES = $(shell mpic++ -showme:compile)
  MPICFLAGS = -DWITHMPI
  MPILFLAGS = $(shell mpic++ -showme:link)
else
  # do nothing..
endif

ifeq ($(USEOMP),true)
  OMPCFLAGS = -fopenmp
  OMPLFLAGS = -fopenmp
else
endif

# Create a build directory for .o object files.
# Crude because this gets run everytime the Makefile
# is parsed. But it works.
$(shell mkdir -p obj)

APPNAME=dvmdostem
LIBDIR=$(SITE_SPECIFIC_LIBS)
INCLUDES=$(SITE_SPECIFIC_INCLUDES) -I external/Fireweed/C++  -I external/Burnup/C++
SOURCES= 	src/TEM.o \
		src/TEMLogger.o \
		src/CalController.o \
		src/ArgHandler.o \
		src/TEMUtilityFunctions.o \
		src/Climate.o \
		src/OutputEstimate.o \
		src/Runner.o \
		src/BgcData.o \
		src/CohortData.o \
		src/EnvData.o \
		src/EnvDataDly.o \
		src/FireData.o \
		src/RestartData.o \
		src/WildFire.o \
		src/DoubleLinkedList.o \
		src/Ground.o \
		src/MineralInfo.o \
		src/Moss.o \
		src/Organic.o \
		src/Snow.o \
		src/SoilParent.o \
		src/Vegetation.o \
		src/CohortLookup.o \
		src/Cohort.o \
		src/Integrator.o \
		src/ModelData.o \
		src/Richards.o \
		src/Snow_Env.o \
		src/Soil_Bgc.o \
		src/Soil_Env.o \
		src/SoilParent_Env.o \
		src/Stefan.o \
		src/TemperatureUpdator.o \
		src/CrankNicholson.o \
		src/tbc-debug-util.o \
		src/Vegetation_Bgc.o \
		src/Vegetation_Env.o \
		src/Layer.o \
		src/MineralLayer.o \
		src/MossLayer.o \
		src/OrganicLayer.o \
		src/ParentLayer.o \
		src/SnowLayer.o \
		src/SoilLayer.o \
		src/FW_Interface.o \
		external/Fireweed/C++/FireweedDeadFuelMoistureFosberg.o \
		external/Fireweed/C++/FireweedFuelModels.o \
		external/Fireweed/C++/FireweedFuelTools.o \
		external/Fireweed/C++/FireweedLiveFuelMoistureGSI.o \
		external/Fireweed/C++/FireweedMessaging.o \
		external/Fireweed/C++/FireweedMetUtils.o \
		external/Fireweed/C++/FireweedRAFireSpread.o \
		external/Fireweed/C++/FireweedStringUtils.o \
		external/Fireweed/C++/FireweedUnits.o \
		external/Fireweed/C++/FireweedUtils.o \
		external/Burnup/C++/BurnupCore.o \
		external/Burnup/C++/BurnupFuelModelInterface.o

OBJECTS =	ArgHandler.o \
		TEMLogger.o \
		CalController.o \
		TEMUtilityFunctions.o \
		Climate.o \
		OutputEstimate.o \
		Runner.o \
		BgcData.o \
		CohortData.o \
		EnvData.o \
		EnvDataDly.o \
		FireData.o \
		RestartData.o \
		WildFire.o \
		DoubleLinkedList.o \
		Ground.o \
		MineralInfo.o \
		Moss.o \
		Organic.o \
		Snow.o \
		SoilParent.o \
		Vegetation.o \
		CohortLookup.o \
		Cohort.o \
		Integrator.o \
		ModelData.o \
		Richards.o \
		Snow_Env.o \
		Soil_Bgc.o \
		Soil_Env.o \
		SoilParent_Env.o \
		Stefan.o \
		CrankNicholson.o \
		tbc-debug-util.o \
		Vegetation_Bgc.o \
		Vegetation_Env.o \
		Layer.o \
		MineralLayer.o \
		MossLayer.o \
		OrganicLayer.o \
		ParentLayer.o \
		SnowLayer.o \
		SoilLayer.o \
		TemperatureUpdator.o \
		FW_Interface.o \
		FireweedDeadFuelMoistureFosberg.o \
		FireweedFuelModels.o \
		FireweedFuelTools.o \
		FireweedLiveFuelMoistureGSI.o \
		FireweedMessaging.o \
		FireweedMetUtils.o \
		FireweedRAFireSpread.o \
		FireweedStringUtils.o \
		FireweedUnits.o \
		FireweedUtils.o \
		BurnupCore.o \
		BurnupFuelModelInterface.o


GIT_SHA := $(shell git describe --abbrev=6 --dirty --always --tags)

TEMOBJ = obj/TEM.o

dvm: $(SOURCES) $(TEMOBJ)
	$(CC) $(SITE_SPECIFIC_LINK_FLAGS) -o $(APPNAME) $(INCLUDES) $(addprefix obj/, $(OBJECTS)) $(TEMOBJ) $(LIBDIR) $(LIBS) $(MPILFLAGS) $(OMPLFLAGS)


lib: $(SOURCES) 
	$(CC) -o libTEM.so -shared $(INCLUDES) $(addprefix obj/, $(OBJECTS)) $(LIBDIR) $(LIBS) $(MPILFLAGS) $(OMPLFLAGS)

CFLAGS += -DGIT_SHA=\"$(GIT_SHA)\"

.cpp.o:
	$(CC) $(CFLAGS) $(MPICFLAGS) $(OMPCFLAGS) $(INCLUDES) $(MPIINCLUDES) $< -o obj/$(notdir $@)

clean:
	rm -f $(OBJECTS) $(APPNAME) TEM.o libTEM.so* *~ obj/*

