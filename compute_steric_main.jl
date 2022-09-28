```
Copyright [2022] JPL

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
```


# --------------------------------------------------------------------------
# Scripts to compute steric sea level and OHC from gridded T and S estimates
# --------------------------------------------------------------------------
printstyled("Loading dependencies...\n",color=:green)
using Distributed
(gethostname() == "MT-315820") & (nprocs() == 1) ? addprocs(10) : nothing
@everywhere using NCDatasets
@everywhere using GibbsSeaWater
@everywhere using SharedArrays
using Dates
using Downloads
using ZipFile
using FTPClient
using DelimitedFiles

printstyled("Loading dependencies done\n",color=:green)
dir_code = homedir()*"/Projects/2022_HOMAGE/"
printstyled("Loading modules...\n",color=:green)
include(dir_code*"mod_clim.jl")
include(dir_code*"mod_compute.jl")
include(dir_code*"mod_dl_cv_BOA.jl")
include(dir_code*"mod_dl_cv_EN4.jl")
include(dir_code*"mod_dl_cv_IAP.jl")
include(dir_code*"mod_dl_cv_SIO.jl")
include(dir_code*"mod_steric_generic.jl")
printstyled("Loading modules done\n",color=:green)
printstyled("  Theads: "*string(Threads.nthreads())*"\n")
printstyled("  Processes: "*string(nprocs())*"\n")

function main()
    settings = def_settings()
    if settings["download_files"] 
        settings["model_name"] == "BOA" ? download_BOA(settings) : nothing
        settings["model_name"] == "EN4_l09" ? download_EN4(settings) : nothing
        settings["model_name"] == "EN4_g10" ? download_EN4(settings) : nothing
        settings["model_name"] == "EN4_c13" ? download_EN4(settings) : nothing
        settings["model_name"] == "EN4_c14" ? download_EN4(settings) : nothing
        settings["model_name"] == "I17" ? printstyled("For I17, use Python script dl_cv_I17.py to download the source data\n",color=:red) : nothing
        settings["model_name"] == "IAP" ? download_IAP(settings) : nothing
        settings["model_name"] == "SIO" ? download_SIO(settings) : nothing
    end

    if settings["convert_files"] 
        settings["model_name"] == "BOA" ? convert_BOA(settings) : nothing
        settings["model_name"] == "EN4_l09" ? convert_EN4(settings) : nothing
        settings["model_name"] == "EN4_g10" ? convert_EN4(settings) : nothing
        settings["model_name"] == "EN4_c13" ? convert_EN4(settings) : nothing
        settings["model_name"] == "EN4_c14" ? convert_EN4(settings) : nothing
        settings["model_name"] == "I17" ? printstyled("For I17, use Python script dl_cv_I17.py to convert the source data\n",color=:red) : nothing
        settings["model_name"] == "IAP" ? convert_IAP(settings) : nothing
        settings["model_name"] == "SIO" ? convert_SIO(settings) : nothing
    end

    if settings["compute_climatology"]
        compute_climatology(settings)
    end
    
    if settings["compute_steric_ohc"]
        compute_steric_ohc(settings)
    end
    return nothing
end

function def_settings()
    # Define all settings
    settings = Dict()
    # Set model name. Currently implemented: 
    # "BOA"         Barnes Optimal 
    # "EN4_l09"     EN4 with fall rate correction from Levitus et al. (2009)  
    # "EN4_g10"     EN4 with fall rate correction from Gouretski et al. (2010)  
    # "EN4_c13"     EN4 with fall rate correction from Cowley et al. (2013)
    # "EN4_c14"     EN4 with fall rate correction from Cheng et al. (2014)  
    # "IAP"         Institute for Atmospheric Physics, Cheng et al. (2016)
    # "I17"         JMA, Ishii et al. (2017)
    # "SIO"         Scripps Institute of Oceanography, Roemmich and Wilson (2009)
    settings["model_name"] = "EN4_c14"

    settings["years"] = [2018:2021...] # Years that will be download, converted and for which steric and OHC will be computed
    settings["years_clim"] = [2018:2020...] # Years that will be used to determine the climatology
    settings["tvec"]  = years2tvec(settings["years"])
    settings["tvec_clim"]  = years2tvec(settings["years_clim"])

    settings["dir_analysis"] = homedir()*"/Data/Steric/"*settings["model_name"]*"/analysis/" # Directory where the files will be downloaded
    settings["dir_analysis_fmt"] = homedir()*"/Data/Steric/"*settings["model_name"]*"/analysis_fmt/" # Directory where the converted monthly T and S files and the climatology will be stored
    settings["fn_savefile"] = homedir()*"/Data/Steric/"*settings["model_name"]*"/postprocessed/"*settings["model_name"]*"_steric_ohc_"*string(settings["years"][1])*"_"*string(settings["years"][end])*".nc" # File where the final results will be saved
    
    # Toggles to set which tasks to perform
    settings["download_files"] = true # Download the raw source data
    settings["convert_files"] = true  # Convert the raw source data to monthly NetCDF files
    settings["compute_climatology"] = true # Compute the climatology
    settings["compute_steric_ohc"] = true # Compute steric and OHC from T and S

    # Do a few checks
    issubset(settings["years_clim"], settings["years"]) ? nothing : printstyled("WARNING: settings['years_clim'] not a subset of settings['years'].\n" ,color=:red, bold=true)
    isdir(settings["dir_analysis"]) ? nothing : printstyled("WARNING: directory "*settings["dir_analysis"]*" does not exist. \n" ,color=:red, bold=true)
    isdir(settings["dir_analysis_fmt"]) ? nothing : printstyled("WARNING: directory "*settings["dir_analysis_fmt"]*" does not exist. \n" ,color=:red, bold=true)
    return settings
end

main()
