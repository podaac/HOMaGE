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

# --------------------------------------------------------------
# Download and convert the Scripps Institute of Oceanography
# temeprature and salinity grids, based on Roemmich & Gilson 
# (2008) into the standard format that can be read by the other 
# scripts
# --------------------------------------------------------------
function download_SIO(settings)
    printstyled(" Downloading SIO data...\n",color=:blue)
    adr_temp = "ftp://kakapo.ucsd.edu/pub/argo/Global_Marine_Argo_Atlas/RG_ArgoClim_Temp.nc"
    adr_psal = "ftp://kakapo.ucsd.edu/pub/argo/Global_Marine_Argo_Atlas/RG_ArgoClim_Psal.nc"
    dl_temp = @task Downloads.download(adr_temp,settings["dir_analysis"]*"RG_ArgoClim_Temp.nc")
    dl_psal = @task Downloads.download(adr_psal,settings["dir_analysis"]*"RG_ArgoClim_Psal.nc")
    schedule(dl_temp)
    schedule(dl_psal)
    wait(dl_temp)
    wait(dl_psal)
    printstyled(" Downloading SIO data done\n",color=:blue)
    return nothing
end

function convert_SIO(settings)
    grid = process_grid_SIO(settings)
    save_grid_info(grid,settings)
    convert_monthly_files_SIO(grid,settings)
end

function process_grid_SIO(settings)
    # Load file and compute all grid variables
    printstyled(" Processing grid properties...\n",color=:blue)
    fh = Dataset(settings["dir_analysis"]*"RG_ArgoClim_Temp.nc","r")
    # Lat,Lon
    θ = fh["LATITUDE"][:]
    ϕ = fh["LONGITUDE"][:]
    # Switch to 0-360 longitude
    shift_lon_idx = [[341:360...];[1:340...]]
    ϕ = mod.(ϕ,360)[shift_lon_idx]

    # Height and pressure
    p = repeat(reshape(fh["PRESSURE"][:],(1,1,:)),length(ϕ),length(θ),1);
    z = @. -convert(Float32,gsw_z_from_p(p, $reshape(θ,(1,:,1))));

    # Top and bottom of each layer
    p_lcl = p[1,1,:]
    plims = zeros(Float32,2,length(p_lcl)) 
    plims[2,1] = (p_lcl[2] - p_lcl[1])/2
    for layer ∈ 2:length(p_lcl)-1
        plims[1,layer] = plims[2,layer-1]
        plims[2,layer] = p_lcl[layer] +(p_lcl[layer+1] - p_lcl[layer])/2
    end
    plims[1,end] = plims[2,end-1]
    plims[2,end] = 2000  

    Δp = repeat(reshape(diff(plims,dims=1)[1,:],(1,1,:)),length(ϕ),length(θ),1);
    ztop = @. gsw_z_from_p($reshape(plims[1,:],(1,1,:)), $reshape(θ,(1,:,1)))
    zbottom = @. gsw_z_from_p($reshape(plims[2,:],(1,1,:)), $reshape(θ,(1,:,1)))
    Δz = repeat((@. ztop - zbottom),length(ϕ),1,1)

    # Area and volume
    A = grid_area(ϕ,θ);
    V = @. A * Δz;

    # Sea-land mask
    slm = (nomissing(fh["BATHYMETRY_MASK"][:],0) .== 1)[shift_lon_idx,:,:]
    @. V[~slm] = 0;
    @. A[~slm[:,:,1]] = 0;
    @. Δz[~slm] = 0;
    @. Δp[~slm] = 0;

    # Store in structure
    grid=grid_struct(ϕ,θ,z,p,Δz,Δp,A,V,slm);
    close(fh)
    return grid
end

function convert_monthly_files_SIO(grid,settings)
    # Convert EN4 files to the generic format
    printstyled(" Converting to generic format...\n",color=:blue)
    shift_lon_idx = [[341:360...];[1:340...]]

    # Read data
    fh_t = Dataset(settings["dir_analysis"]*"RG_ArgoClim_Temp.nc","r")
    fh_s = Dataset(settings["dir_analysis"]*"RG_ArgoClim_Psal.nc","r")

    traw = fh_t["TIME"][:];
    sio_tvec = zeros(Int32,length(traw),2)
    sio_tvec[:,1] = fld.(traw,12) .+ 2004
    sio_tvec[:,2] = mod.(traw,12) .+ 0.5

    temp_mean = nomissing(fh_t["ARGO_TEMPERATURE_MEAN"][:],0)[shift_lon_idx,:,:];
    psal_mean = nomissing(fh_s["ARGO_SALINITY_MEAN"][:],0)[shift_lon_idx,:,:];
    close(fh_t)
    close(fh_s)
    
    @sync @distributed for tstep in axes(settings["tvec"],1)
        yr = settings["tvec"][tstep,1]
        mnth = settings["tvec"][tstep,2]
        println("  Converting year "*string(yr)*" month "*string(mnth)*"...")
        read_idx = findfirst((sio_tvec[:,1] .== yr) .& (sio_tvec[:,2] .== mnth))
        fh_t = Dataset(settings["dir_analysis"]*"RG_ArgoClim_Temp.nc","r")
        fh_s = Dataset(settings["dir_analysis"]*"RG_ArgoClim_Psal.nc","r")
        CT = nomissing(fh_t["ARGO_TEMPERATURE_ANOMALY"][:,:,:,read_idx],0)[shift_lon_idx,:,:];
        SA = nomissing(fh_s["ARGO_SALINITY_ANOMALY"][:,:,:,read_idx],0)[shift_lon_idx,:,:];
        close(fh_t)
        close(fh_s)
     
        # Convert to SA,CT,ρ
        @. SA = grid.slm * gsw_sa_from_sp(SA+psal_mean,grid.p,$reshape(grid.ϕ,(:,1,1)),$reshape(grid.θ,(1,:,1)));
        @. CT = grid.slm * gsw_ct_from_t(SA,CT+temp_mean,grid.p);  
        ρ = @. grid.slm * convert(Float32,gsw_rho(SA,CT,grid.p));
        data_monthly = data_struct(SA,CT,ρ)
        save_analysis_fmt(settings["dir_analysis_fmt"],settings["model_name"],yr,mnth,grid,data_monthly) # Save file
    end
    return nothing
end