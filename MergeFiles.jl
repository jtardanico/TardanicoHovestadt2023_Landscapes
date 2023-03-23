# Data file merge script

# For use combining simulation replicate output files.into a single data file
# Group files by simulation parameters.

println("Loading packages")

using DataFrames
using CSV
#using Parsers
using ArgParse

#------------------------------------
# Function definitions
#------------------------------------

# Read in read_arguments from the shell script
function read_args()
    s=ArgParseSettings()
    @add_arg_table s begin
        "--parasource","-n"
            help = "parameter dictionary source file"
    end
    return parse_args(s)
end

# Outputs diagnostic info to console
function df_info(data::DataFrame)
    println("Data file dimensions")
    println("---------------------")
    println("cols: ",length(data[1,1:end]))
    println("rows: ",length(data[1:end,1]))
    println("")
    println("Missing values and NaNs")
    println("-----------------------")
    for col in names(data)
        println("Column $col")
        println("Length:")
        println(length(data[1:end,Symbol(col)]))
        println("Missing:")
        println(length(filter(ismissing,data[:,Symbol(col)])))
        if length(filter(x->x==true,isa.(data[:,Symbol(col)],Number))) > 0
            println("Nans:")
            println(length(filter(isnan,data[:,Symbol(col)])))
        end
    end
    println("")
end

# Gets the replicate number of a data file
function extract_repli_num(filename::String,regex::Regex)
    m = match(regex,filename)
    y = m.match
    y = replace(y,"_"=>"")
    y = replace(y,".txt"=>"")
    y = replace(y,"trend"=>"")
    y = replace(y,"patches"=>"")
    println("Replicate $y")
    y = parse(Int,y)
    return y
end

#----------------------------------
# Main Script
#----------------------------------

println("Starting main script")
args = read_args()
println("Parasource = ",args["parasource"])
par = include(args["parasource"])
dir = par["dir"]
println("File directory = ",dir)

#---------------------------------------
# Use if running directly from IDE
#---------------------------------------

#using DataFrames
#using CSV
#using ArgParse

#dir = ""


#---------------------------------------

infiles = readdir(dir)

filekeys = copy(infiles)

s = r"_\d+\.txt"
s2 = r"_\d+trend\.txt"
s3 = r"_\d+patches\.txt"


println("Getting file name keywords")
for i in 1:length(filekeys)
    #filekeys[i] = replace(filekeys[i],dir=>"")
    filekeys[i] = replace(filekeys[i],s=>"")
    filekeys[i] = replace(filekeys[i],s2=>"")
    filekeys[i] = replace(filekeys[i],s3=>"")
end

unique!(filekeys)
println(filekeys)

filekeys = filekeys[occursin.(Regex("console_output"),filekeys).==false]
#filekeys = filekeys[occursin.(Regex("merged"),filekeys).==false]

println("Merging data files")
println("---------------------")
for i in 1:length(filekeys)
    # Merge whole landscape data files
    println(filekeys[i])
    println("Individual data")
    infiles2 = infiles[occursin.(Regex("trend"),infiles).==false]
    infiles2 = infiles2[occursin.(Regex("patches"),infiles2).==false]
    infiles2 = infiles2[occursin.(Regex("console_output"),infiles2).==false]
    infiles2 = infiles2[occursin.(Regex(filekeys[i]),infiles2).==true]
    outfile = string(dir,"merged",filekeys[i],".txt")
    println(length(infiles2))
    #println(outfile)
    for j in 1:length(infiles2)
        println("file $j")
        println(infiles2[j])
        file = DataFrame(CSV.File(string(dir,infiles2[j])))
        #df_info(file)
        file.Replicate .= extract_repli_num(infiles2[j],s)

        if j==1

            CSV.write(outfile,file,append=false)
        else
            CSV.write(outfile,file,append=true)
        end
    end
    # Merge trend data files
    #println(filekeys[i]," trend")
    println("Trend data")
    infiles2 = infiles[occursin.(Regex("trend"),infiles).==true]
    infiles2 = infiles2[occursin.(Regex("patches"),infiles2).==false]
    infiles2 = infiles2[occursin.(Regex(filekeys[i]),infiles2).==true]
    outfile = string(dir,"merged",filekeys[i],"_trend.txt")
    println(length(infiles2))
    println(infiles2)
    #println(outfile)
    for j in 1:length(infiles2)
        println("file $j")
        println(infiles2[j])
        file = DataFrame(CSV.File(string(dir,infiles2[j])))
        #df_info(file)
        file.Replicate .= extract_repli_num(infiles2[j],s2)
        if j==1
            CSV.write(outfile,file,append=false)
        else
            CSV.write(outfile,file,append=true)
        end
    end
    #println(filekeys[i]," patches")
    println("Patch data")
    infiles2 = infiles[occursin.(Regex("trend"),infiles).==false]
    infiles2 = infiles2[occursin.(Regex("patches"),infiles2).==true]
    infiles2 = infiles2[occursin.(Regex(filekeys[i]),infiles2).==true]
    outfile = string(dir,"merged",filekeys[i],"_patches.txt")
    println(length(infiles2))
    println(infiles2)
    for j in 1:length(infiles2)
        println("file $j")
        println(infiles2[j])
        file = DataFrame(CSV.File(string(dir,infiles2[j])))
        #df_info(file)
        file.Replicate .= extract_repli_num(infiles2[j],s3)
        if j==1
            CSV.write(outfile,file,append=false)
        else
            CSV.write(outfile,file,append=true)
        end
    end
    println("end script")
end
