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

# -----------------------------------
# Compute a climatology for SA, CT, ρ
# At this point, it only incorporates 
# a single value, so no monthly
# climatology
# -----------------------------------
function compute_climatology(settings)
    printstyled(" Computing climatology...\n",color=:blue)
    grid = load_grid_info(settings); # Read grid file

    printstyled("  Salinity...\n",color=:blue)
    SA_c = @distributed (+) for tstep in axes(settings["tvec_clim"],1)
    read_field(settings["dir_analysis_fmt"] * settings["model_name"]*"_"*string(settings["tvec_clim"][tstep,1])*"_"*lpad(settings["tvec_clim"][tstep,2],2,"0")*".nc","SA");
    end;

    printstyled("  Temperature...\n",color=:blue)
    CT_c = @distributed (+) for tstep in axes(settings["tvec_clim"],1)
    read_field(settings["dir_analysis_fmt"] * settings["model_name"]*"_"*string(settings["tvec_clim"][tstep,1])*"_"*lpad(settings["tvec_clim"][tstep,2],2,"0")*".nc","CT");
    end;

    printstyled("  Density...\n",color=:blue)
    ρ_c = @distributed (+) for tstep in axes(settings["tvec_clim"],1)
    read_field(settings["dir_analysis_fmt"] * settings["model_name"]*"_"*string(settings["tvec_clim"][tstep,1])*"_"*lpad(settings["tvec_clim"][tstep,2],2,"0")*".nc","rho");
    end;

    # Divide by the number of time steps to get the mean instead of the sum 
    SA_c ./= size(settings["tvec_clim"],1);
    CT_c ./= size(settings["tvec_clim"],1);
    ρ_c ./= size(settings["tvec_clim"],1);

    printstyled("  Saving...\n",color=:blue)
    data_clim = data_struct(SA_c,CT_c,ρ_c)
    save_climatology(data_clim,grid,settings)
    printstyled(" Computing climatology done\n",color=:blue)
    return nothing
end

@everywhere function read_field(fn,varname)
    fh = Dataset(fn,"r")
    field = fh[varname][:] :: Array{Float32,3};
    close(fh)
    return field
end

function save_climatology(data_clim,grid,settings)
    fn = settings["dir_analysis_fmt"] * settings["model_name"]*"_climatology.nc"
    fh = Dataset(fn,"c")
    defDim(fh,"lon",length(grid.ϕ))
    defDim(fh,"lat",length(grid.θ))
    defDim(fh,"lvl",size(grid.z,3))
    defVar(fh,"SA",data_clim.SA,("lon","lat","lvl"),deflatelevel=1)
    defVar(fh,"CT",data_clim.CT,("lon","lat","lvl",),deflatelevel=1)
    defVar(fh,"rho",data_clim.ρ,("lon","lat","lvl"),deflatelevel=1)
    close(fh)
    return nothing
end

@everywhere function load_climatology(settings)
    fh = Dataset(settings["dir_analysis_fmt"] * settings["model_name"]*"_climatology.nc","r")
    SA = fh["SA"][:] :: Array{Float32,3};
    CT = fh["CT"][:] :: Array{Float32,3};
    ρ = fh["rho"][:] :: Array{Float32,3}
    data_clim = data_struct(SA,CT,ρ);
    close(fh)
    return(data_clim)
end