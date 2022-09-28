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

### General steric functions and structures
# ----------
# Structures
# ----------
@everywhere struct grid_struct
    # Structure to hold all grid variables
    ϕ :: Array{Float32,1} # Longitude (°)
    θ :: Array{Float32,1} # Latitude (°)
    z :: Array{Float32,3} # Depth (m)
    p :: Array{Float32,3} # Pressure (Pa)
    Δz :: Array{Float32,3} # Height of each grid cell (m)
    Δp :: Array{Float32,3} # Pressure difference of each grid cell (Pa)
    A :: Array{Float32,2} # Area of each grid cell (m²)
    V :: Array{Float32,3} # Volume of each grid cell (m³)
    slm :: Array{Bool,3} # Sea-land mask
end

@everywhere struct data_struct
    # Structure to hold all data variables
    SA :: Array{Float32,3} # Pressure difference of each grid cell (Pa)
    CT :: Array{Float32,3} # Area of each grid cell (m²)
    ρ :: Array{Float32,3} # Volume of each grid cell (m³)
end

# ------------------
# Grid file functions
# ------------------
function save_grid_info(grid,settings)
    fn = settings["dir_analysis_fmt"] * settings["model_name"]*"_grid.nc"
    fh = Dataset(fn,"c")
    defDim(fh,"lon",length(grid.ϕ))
    defDim(fh,"lat",length(grid.θ))
    defDim(fh,"lvl",size(grid.z,3))
    defVar(fh,"lvl",[1:size(grid.z,3)...],("lvl",),deflatelevel=5)
    defVar(fh,"lon",grid.ϕ,("lon",),deflatelevel=5)
    defVar(fh,"lat",grid.θ,("lat",),deflatelevel=5)

    defVar(fh,"z",grid.z,("lon","lat","lvl"),deflatelevel=5)
    defVar(fh,"p",grid.p,("lon","lat","lvl",),deflatelevel=5)

    defVar(fh,"Dz",grid.Δz,("lon","lat","lvl"),deflatelevel=5)
    defVar(fh,"Dp",grid.Δp,("lon","lat","lvl",),deflatelevel=5)

    defVar(fh,"A",grid.A,("lon","lat"),deflatelevel=5)
    defVar(fh,"V",grid.V,("lon","lat","lvl",),deflatelevel=5)
    defVar(fh,"slm",convert.(Int8,grid.slm),("lon","lat","lvl",),deflatelevel=5)
    close(fh)
    return nothing
end

@everywhere function load_grid_info(settings)
    fn = settings["dir_analysis_fmt"] * settings["model_name"]*"_grid.nc"
    fh = Dataset(fn,"r")
    ϕ = fh["lon"][:]
    θ = fh["lat"][:]
    z = fh["z"][:]
    p = fh["p"][:]
    Δz = fh["Dz"][:]
    Δp = fh["Dp"][:]
    A = fh["A"][:]
    V = fh["V"][:]
    slm = convert.(Bool, fh["slm"][:])
    close(fh)
    grid=grid_struct(ϕ,θ,z,p,Δz,Δp,A,V,slm);
    return grid
end

@everywhere function grid_area(ϕ,θ)
    gridsize = abs(θ[2]-θ[1])
    area = @. deg2rad(gridsize) * (sind(θ+gridsize/2)-sind(θ-gridsize/2)) * 6371000^2
    return repeat(area',length(ϕ))
end

# -------------------
# Data file functions
# -------------------
@everywhere function save_analysis_fmt(dir,model_name,yr,mnth,grid,data_monthly)
    fh = Dataset(dir * model_name*"_"*string(yr)*"_"*lpad(mnth,2,"0")*".nc","c")
    defDim(fh,"lon",length(grid.ϕ))
    defDim(fh,"lat",length(grid.θ))
    defDim(fh,"lvl",size(grid.z,3))
    defVar(fh,"lvl",[1:size(grid.z,3)...],("lvl",),deflatelevel=1)
    defVar(fh,"lon",grid.ϕ,("lon",),deflatelevel=1)
    defVar(fh,"lat",grid.θ,("lat",),deflatelevel=1)
    defVar(fh,"SA",data_monthly.SA,("lon","lat","lvl"),deflatelevel=1)
    defVar(fh,"CT",data_monthly.CT,("lon","lat","lvl",),deflatelevel=1)
    defVar(fh,"rho",data_monthly.ρ,("lon","lat","lvl"),deflatelevel=1)
    close(fh)
    return nothing
end

@everywhere function read_analysis_fmt(dir,model_name,yr,mnth)
    fh = Dataset(dir * model_name*"_"*string(yr)*"_"*lpad(mnth,2,"0")*".nc","r")
    SA = fh["SA"][:] :: Array{Float32,3}
    CT = fh["CT"][:] :: Array{Float32,3}
    ρ = fh["rho"][:] :: Array{Float32,3}
    data_monthly = data_struct(SA,CT,ρ);
    close(fh)
    return(data_monthly)
end

function years2tvec(years)
    tvec  = zeros(Int32,12*length(years),2)
    for (idx,yr) in enumerate(years)
        @. tvec[12*(idx-1)+1:12*idx,1] = yr
        @. tvec[12*(idx-1)+1:12*idx,2] = [1:12...]
    end
    return tvec
end