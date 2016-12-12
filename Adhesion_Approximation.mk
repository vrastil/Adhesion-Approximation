##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=Adhesion_Approximation
ConfigurationName      :=Debug
WorkspacePath          := "/home/vrastil/Dropbox/Argonne/Adhesion_Approximation"
ProjectPath            := "/home/vrastil/Dropbox/Argonne/Adhesion_Approximation"
IntermediateDirectory  :=./Debug
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=Michal Vrastil
Date                   :=06/12/16
CodeLitePath           :="/home/vrastil/.codelite"
LinkerName             :=/usr/bin/x86_64-linux-gnu-g++
SharedObjectLinkerName :=/usr/bin/x86_64-linux-gnu-g++ -shared -fPIC
ObjectSuffix           :=.o
DependSuffix           :=.o.d
PreprocessSuffix       :=.i
DebugSwitch            :=-g 
IncludeSwitch          :=-I
LibrarySwitch          :=-l
OutputSwitch           :=-o 
LibraryPathSwitch      :=-L
PreprocessorSwitch     :=-D
SourceSwitch           :=-c 
OutputFile             :=$(IntermediateDirectory)/$(ProjectName)
Preprocessors          :=
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E
ObjectsFileList        :="Adhesion_Approximation.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  -O0
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). 
IncludePCH             := 
RcIncludePath          := 
Libs                   := 
ArLibs                 :=  
LibPath                := $(LibraryPathSwitch). $(LibraryPathSwitch). $(LibraryPathSwitch)Debug 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := /usr/bin/x86_64-linux-gnu-ar rcu
CXX      := /usr/bin/x86_64-linux-gnu-g++
CC       := /usr/bin/x86_64-linux-gnu-gcc
CXXFLAGS :=  -g -Wall $(Preprocessors)
CFLAGS   :=   $(Preprocessors)
ASFLAGS  := 
AS       := /usr/bin/x86_64-linux-gnu-as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=$(IntermediateDirectory)/src_output.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_mod_frozen_pot.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_approximations.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_cmd_line.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_mesh.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_main.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_grid_fce.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_CBRNG_Random.cpp$(ObjectSuffix) $(IntermediateDirectory)/src_random.cpp$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: $(OutputFile)

$(OutputFile): $(IntermediateDirectory)/.d $(Objects) 
	@$(MakeDirCommand) $(@D)
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(LinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)

MakeIntermediateDirs:
	@test -d ./Debug || $(MakeDirCommand) ./Debug


$(IntermediateDirectory)/.d:
	@test -d ./Debug || $(MakeDirCommand) ./Debug

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/src_output.cpp$(ObjectSuffix): src/output.cpp $(IntermediateDirectory)/src_output.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/vrastil/Dropbox/Argonne/Adhesion_Approximation/src/output.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_output.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_output.cpp$(DependSuffix): src/output.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_output.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_output.cpp$(DependSuffix) -MM "src/output.cpp"

$(IntermediateDirectory)/src_output.cpp$(PreprocessSuffix): src/output.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_output.cpp$(PreprocessSuffix) "src/output.cpp"

$(IntermediateDirectory)/src_mod_frozen_pot.cpp$(ObjectSuffix): src/mod_frozen_pot.cpp $(IntermediateDirectory)/src_mod_frozen_pot.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/vrastil/Dropbox/Argonne/Adhesion_Approximation/src/mod_frozen_pot.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_mod_frozen_pot.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_mod_frozen_pot.cpp$(DependSuffix): src/mod_frozen_pot.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_mod_frozen_pot.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_mod_frozen_pot.cpp$(DependSuffix) -MM "src/mod_frozen_pot.cpp"

$(IntermediateDirectory)/src_mod_frozen_pot.cpp$(PreprocessSuffix): src/mod_frozen_pot.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_mod_frozen_pot.cpp$(PreprocessSuffix) "src/mod_frozen_pot.cpp"

$(IntermediateDirectory)/src_approximations.cpp$(ObjectSuffix): src/approximations.cpp $(IntermediateDirectory)/src_approximations.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/vrastil/Dropbox/Argonne/Adhesion_Approximation/src/approximations.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_approximations.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_approximations.cpp$(DependSuffix): src/approximations.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_approximations.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_approximations.cpp$(DependSuffix) -MM "src/approximations.cpp"

$(IntermediateDirectory)/src_approximations.cpp$(PreprocessSuffix): src/approximations.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_approximations.cpp$(PreprocessSuffix) "src/approximations.cpp"

$(IntermediateDirectory)/src_cmd_line.cpp$(ObjectSuffix): src/cmd_line.cpp $(IntermediateDirectory)/src_cmd_line.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/vrastil/Dropbox/Argonne/Adhesion_Approximation/src/cmd_line.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_cmd_line.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_cmd_line.cpp$(DependSuffix): src/cmd_line.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_cmd_line.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_cmd_line.cpp$(DependSuffix) -MM "src/cmd_line.cpp"

$(IntermediateDirectory)/src_cmd_line.cpp$(PreprocessSuffix): src/cmd_line.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_cmd_line.cpp$(PreprocessSuffix) "src/cmd_line.cpp"

$(IntermediateDirectory)/src_mesh.cpp$(ObjectSuffix): src/mesh.cpp $(IntermediateDirectory)/src_mesh.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/vrastil/Dropbox/Argonne/Adhesion_Approximation/src/mesh.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_mesh.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_mesh.cpp$(DependSuffix): src/mesh.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_mesh.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_mesh.cpp$(DependSuffix) -MM "src/mesh.cpp"

$(IntermediateDirectory)/src_mesh.cpp$(PreprocessSuffix): src/mesh.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_mesh.cpp$(PreprocessSuffix) "src/mesh.cpp"

$(IntermediateDirectory)/src_main.cpp$(ObjectSuffix): src/main.cpp $(IntermediateDirectory)/src_main.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/vrastil/Dropbox/Argonne/Adhesion_Approximation/src/main.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_main.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_main.cpp$(DependSuffix): src/main.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_main.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_main.cpp$(DependSuffix) -MM "src/main.cpp"

$(IntermediateDirectory)/src_main.cpp$(PreprocessSuffix): src/main.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_main.cpp$(PreprocessSuffix) "src/main.cpp"

$(IntermediateDirectory)/src_grid_fce.cpp$(ObjectSuffix): src/grid_fce.cpp $(IntermediateDirectory)/src_grid_fce.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/vrastil/Dropbox/Argonne/Adhesion_Approximation/src/grid_fce.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_grid_fce.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_grid_fce.cpp$(DependSuffix): src/grid_fce.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_grid_fce.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_grid_fce.cpp$(DependSuffix) -MM "src/grid_fce.cpp"

$(IntermediateDirectory)/src_grid_fce.cpp$(PreprocessSuffix): src/grid_fce.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_grid_fce.cpp$(PreprocessSuffix) "src/grid_fce.cpp"

$(IntermediateDirectory)/src_CBRNG_Random.cpp$(ObjectSuffix): src/CBRNG_Random.cpp $(IntermediateDirectory)/src_CBRNG_Random.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/vrastil/Dropbox/Argonne/Adhesion_Approximation/src/CBRNG_Random.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_CBRNG_Random.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_CBRNG_Random.cpp$(DependSuffix): src/CBRNG_Random.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_CBRNG_Random.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_CBRNG_Random.cpp$(DependSuffix) -MM "src/CBRNG_Random.cpp"

$(IntermediateDirectory)/src_CBRNG_Random.cpp$(PreprocessSuffix): src/CBRNG_Random.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_CBRNG_Random.cpp$(PreprocessSuffix) "src/CBRNG_Random.cpp"

$(IntermediateDirectory)/src_random.cpp$(ObjectSuffix): src/random.cpp $(IntermediateDirectory)/src_random.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/vrastil/Dropbox/Argonne/Adhesion_Approximation/src/random.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_random.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_random.cpp$(DependSuffix): src/random.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_random.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/src_random.cpp$(DependSuffix) -MM "src/random.cpp"

$(IntermediateDirectory)/src_random.cpp$(PreprocessSuffix): src/random.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_random.cpp$(PreprocessSuffix) "src/random.cpp"


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Debug/


