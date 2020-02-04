# Master IPM sourcing script

# See each sub-script for notes on specifics

source('Ana_Israel_IPM/R/01_Utility_Functions_and_Dependencies.R')
sourceCpp('Ana_Israel_IPM/Cpp/01_Utility_funs_cpp.cpp')
source('Ana_Israel_IPM/R/02_Read_Data.R')
source('Ana_Israel_IPM/R/03_VR_Model_Selection.R')
source('Ana_Israel_IPM/R/04_IPM_Implementation.R')
source('Ana_Israel_IPM/R/05_Sens_and_Elas.R')
source('Ana_Israel_IPM/R/06_Bootstrap.R')
source('Ana_Israel_IPM/R/07_Figures.R')
