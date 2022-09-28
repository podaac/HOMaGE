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

# Compute steric sea level and OHC from SA, CT, ρ anomalies
function compute_steric_ohc(settings)
    printstyled(" Computing steric and ocean heat content...\n",color=:blue)
    printstyled("  Reading grid info and climatology...\n",color=:blue)
    grid      = load_grid_info(settings);   # Load grid file
    data_clim = load_climatology(settings); # Load climatology

    # These are the 3D (ϕ,θ,time) fields where we store gridded:
    printstyled("  Allocating arrays...\n",color=:blue)
    η    = SharedArray{Float32,3}(length(grid.ϕ),length(grid.θ),size(settings["tvec"],1));    # steric (η) (m)
    η_SA = SharedArray{Float32,3}(length(grid.ϕ),length(grid.θ),size(settings["tvec"],1));    # halosteric (η_SA) (m)
    η_CT = SharedArray{Float32,3}(length(grid.ϕ),length(grid.θ),size(settings["tvec"],1));    # thermosteric (η_CT) (m)
    ζ    = SharedArray{Float32,3}(length(grid.ϕ),length(grid.θ),size(settings["tvec"],1));    # Ocean heat content (ζ) (J/grid cell) (We might want to change this to J/m^2 ???)

    # Define two functions that compute ohc and steric anomalies for each lat/lon/depth point 
    @everywhere ohc_fun = (V,CT,ρ,CT_clim,ρ_clim) -> V * 3991.8679571196f0 * ((CT*ρ) - (CT_clim*ρ_clim))
    @everywhere str_fun = (Δz,ρ,ρ_clim) -> Δz * (ρ_clim/(ρ+1.0f-10)-1)

    # Compute gridded steric/OHC estimates
    printstyled("  Computing gridded steric and OHC...\n",color=:blue)
    @sync @distributed for tstep in axes(settings["tvec"],1)
        yr   = settings["tvec"][tstep,1]
        mnth = settings["tvec"][tstep,2]
        println("  Processing year "*string(yr)*" month "*string(mnth)*"...")
        data_monthly = read_analysis_fmt(settings["dir_analysis_fmt"],settings["model_name"],yr,mnth); # Read data

        η[:,:,tstep] = mapreduce(str_fun,+,grid.Δz,data_monthly.ρ,data_clim.ρ,dims=3);
        ρ = @. grid.slm * convert(Float32,gsw_rho(data_monthly.SA,data_clim.CT,grid.p));
        η_SA[:,:,tstep] = mapreduce(str_fun,+,grid.Δz,ρ,data_clim.ρ,dims=3);        
        @. ρ =  grid.slm * convert(Float32,gsw_rho(data_clim.SA,data_monthly.CT,grid.p));
        η_CT[:,:,tstep] = mapreduce(str_fun,+,grid.Δz,ρ,data_clim.ρ,dims=3);
        ζ[:,:,tstep] = mapreduce(ohc_fun,+,grid.V,data_monthly.CT,data_monthly.ρ,data_clim.CT,data_clim.ρ,dims=3)[:,:,1];
    end

    printstyled("  Computing global-mean values...\n",color=:blue)
    A_ocn  = sum(@. grid.A * grid.slm[:,:,1]); # Ocean area
    ζ_g    = sum(ζ,dims=(1,2))[1,1,:];
    η_g    = sum((@. η * grid.A),dims=(1,2))[1,1,:] ./ A_ocn;
    η_CT_g = sum((@. η_CT * grid.A),dims=(1,2))[1,1,:] ./ A_ocn;
    η_SA_g = sum((@. η_SA * grid.A),dims=(1,2))[1,1,:] ./ A_ocn;

    @. η -= $reshape(η_SA_g,(1,1,:));    # Remove global-mean halosteric from gridded steric due to salinity drift

    printstyled("  Saving results...\n",color=:blue)
    save_steric_ohc(η,η_SA,η_CT,ζ,ζ_g, η_CT_g, η_SA_g, η_g,grid,settings)      # Save everything
    printstyled(" Computing steric and ocean heat content done\n",color=:blue)
    return nothing
end


function save_steric_ohc(η,η_SA,η_CT,ζ,ζ_g, η_CT_g, η_SA_g, η_g,grid,settings)
    # Save the data
    fh = Dataset(settings["fn_savefile"],"c")
    fh.attrib["data_source"] = settings["model_name"]
    tval = @. DateTime(settings["tvec"][:,1],settings["tvec"][:,2],15)
    defDim(fh,"time",length(tval))
    defDim(fh,"lon",length(grid.ϕ))
    defDim(fh,"lat",length(grid.θ))
    defDim(fh,"depth",size(grid.z,3))
    defVar(fh,"time",tval,("time",),deflatelevel=5)
    defVar(fh,"lon",grid.ϕ,("lon",),fillvalue=-9999,deflatelevel=5,attrib = Dict("long_name" => "Longitude","units" => "Degrees"))
    defVar(fh,"lat",grid.θ,("lat",),fillvalue=-9999,deflatelevel=5,attrib = Dict("long_name" => "Latitude","units" => "Degrees"))
    defVar(fh,"slm",Int8,("lon","lat",),deflatelevel=5,attrib = Dict("long_name" => "Sea-land mask","description" => "1 = Sea, 0 = land"))[:] = grid.slm[:,:,1];
    defVar(fh,"ohc_2d",ζ,("lon","lat","time"),deflatelevel=5,attrib = Dict("long_name" => "Ocean heat content (gridded)","units" => "Joules"))
    defVar(fh,"thermosteric_2d",η_CT,("lon","lat","time"),deflatelevel=5,attrib = Dict("long_name" => "Thermosteric sea-level anomaly (gridded)","units" => "m"))
    defVar(fh,"halosteric_2d",η_SA,("lon","lat","time"),deflatelevel=5,attrib = Dict("long_name" => "Halosteric sea-level anomaly (gridded)","units" => "m"))
    defVar(fh,"totalsteric_2d",η,("lon","lat","time"),deflatelevel=5,attrib = Dict("long_name" => "Steric sea-level anomaly (gridded)","units" => "m"))
    defVar(fh,"ohc_ts",ζ_g,("time",),deflatelevel=5,attrib = Dict("long_name" => "Global ocean heat content","units" => "J"))
    defVar(fh,"thermosteric_ts",η_CT_g,("time",),deflatelevel=5,attrib = Dict("long_name" => "Global-mean thermosteric sea-level anomaly","units" => "m"))
    defVar(fh,"halosteric_ts",η_SA_g,("time",),deflatelevel=5,attrib = Dict("long_name" => "Global-mean halosteric sea-level anomaly","units" => "m"))
    defVar(fh,"totalsteric_ts",η_g,("time",),deflatelevel=5,attrib = Dict("long_name" => "Global-mean steric sea-level anomaly","units" => "m"))
    close(fh)
    return nothing
end
