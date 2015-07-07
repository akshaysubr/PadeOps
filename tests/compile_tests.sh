# Compile test_cd10.F90
ifort -warn all ../src/utilities/kind_parameters.F90 ../src/utilities/constants.F90 ../src/utilities/cd10.F90 test_cd10.F90 -o test_cd10

# Compile test_cd06.F90
ifort -warn all ../src/utilities/kind_parameters.F90 ../src/utilities/constants.F90 ../src/utilities/cd06.F90 test_cd06.F90 -o test_cd06

# Compile test_ffts.F90
ifort -warn all -mkl ../src/utilities/kind_parameters.F90 ../src/utilities/constants.F90 ../src/utilities/ffts.F90 test_ffts.F90 -o test_ffts

# Compile test_dcts.F90
ifort -warn all -mkl ../src/utilities/kind_parameters.F90 ../src/utilities/constants.F90 ../src/utilities/dcts.F90 test_dcts.F90 -o test_dcts
