# Template Shell Script

cd /home/ubuntu/julia-1.1.1/bin/

echo "Starting script..."

# Main program files
source=/home/ubuntu/sims/model1/LandscapesModel1.jl # The main simulation program
para=/home/ubuntu/sims/model1/ParaDict_Int_Test.jl # Parameter file
merge=/home/ubuntu/sims/model1/MergeFiles.jl # Script for merging output files

# Landscape files: Temperature
temp1=/home/ubuntu/landscapefiles/testing/test_AC_temp_20_0.5_0.5_1.txt
temp2=/home/ubuntu/landscapefiles/testing/test_AC_temp_20_0.5_0.5_2.txt
temp3=/home/ubuntu/landscapefiles/testing/test_AC_temp_20_0.5_0.5_3.txt

# Landscape files: "Environment"
env1=/home/ubuntu/landscapefiles/testing/test_AC_env_20_0.5_0.5_1.txt
env2=/home/ubuntu/landscapefiles/testing/test_AC_env_20_0.5_0.5_2.txt
env3=/home/ubuntu/landscapefiles/testing/test_AC_env_20_0.5_0.5_3.txt

# Landscape files: Uniform landscape
unif=/home/ubuntu/landscapefiles/testing/testing_uniformenv.txt

echo "source = $source"
echo "para = $para"

# Simulation runs. See comments in read_arguments() in LandscapeModel1.jl for arguement key.

echo "Uniform env"

./julia -p 1 $source -n $para -t $temp1  -e $unif -s 1 -b false -a true -c 0.5 -i false -o true #1
./julia -p 1 $source -n $para -t $temp2  -e $unif -s 1 -b false -a true -c 0.5 -i false -o true #2
./julia -p 1 $source -n $para -t $temp3  -e $unif -s 1 -b false -a true -c 0.5 -i false -o true #3

echo "Autocorrelated env"

./julia -p 1 $source -n $para -t $temp1  -e $env1 -s 1 -b false -a true -c 0.5 -i true -d 0.5 #11
./julia -p 1 $source -n $para -t $temp2  -e $env2 -s 1 -b false -a true -c 0.5 -i true -d 0.5 #12
./julia -p 1 $source -n $para -t $temp3  -e $env3 -s 1 -b false -a true -c 0.5 -i true -d 0.5 #13

echo "Merging output files"

./julia -p 1 $merge -n $para

echo "Script complete."
