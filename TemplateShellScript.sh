# Tlate Shell Script

cd /julia_directory/julia-1.1.1/bin/ # Location where you installed Julia 1.1.1

echo "Starting script..."

# Main program files
source=/program_directory/LandscapesModel1.jl # This main simulation program file
para=/parameter_config_directory/parameter_config_file.jl # The parameter configuration file
merge=/program_directory/MergeFiles.jl # Script for merging output files

# Landscape files for T attribute
T1=/my_landscape_file_directory/landscapefile1.txt # The directory where your landscape files are located
T2=/my_landscape_file_directory/landscapefile2.txt
T3=/my_landscape_file_directory/landscapefile3.txt

# Landscape files for H attribute
H1=/my_landscape_file_directory/landscapefile1.txt
H2=/my_landscape_file_directory/landscapefile2.txt
H3=/my_landscape_file_directory/landscapefile3.txt


echo "source = $source"
echo "para = $para"

# Simulation runs. See comments in read_arguments() in LandscapeModel1.jl for arguement key.

echo "Running simulations"

./julia -p 1 $source -n $para -t $T1  -e $H -s 1 -b false -a true -c 0.5 -i false -o true # Simulation run #1
./julia -p 1 $source -n $para -t $T2  -e $H -s 1 -b false -a true -c 0.5 -i false -o true # Simulation run #1
./julia -p 1 $source -n $para -t $T3  -e $H -s 1 -b false -a true -c 0.5 -i false -o true # Simulation run #1


echo "Merging output files"

./julia -p 1 $merge -n $para # Combine the output files together (Optional step)

echo "Script complete."
