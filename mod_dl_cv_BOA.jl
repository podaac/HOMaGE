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

# ----------------------------------------------------------------------
# Download and convert the CSIO / BOA temperature and salinity 
# grids, based on Li et al. (2017) into the standard format that 
# can be read by the other scripts.
#
# Li, H., Xu, F., Zhou, W., Wang, D., Wright, J. S., Liu, Z., & Lin, Y. 
# (2017). Development of a global gridded Argo data set with Barnes 
# successive corrections. Journal of Geophysical Research: Oceans,
# 122(2), 866–889. https://doi.org/10.1002/2016JC012285
#
# Data is downloaded from Chinese server. Might not work while on lab
# network !!!
# ----------------------------------------------------------------------
function download_BOA(settings)
    ftp = FTP(hostname="data.argo.org.cn/",implicit=true)
    cd(ftp,"/pub/ARGO/BOA_Argo/NetCDF/")
    flist = readdir(ftp)
    close(ftp)
    tlist = []
    for f in flist
        if length(f) == 19 && all(isnumeric(c) for c in f[10:13]) && all(isnumeric(c) for c in f[15:16]) && (parse(Int,f[10:13]) in settings["tvec"][:,1])
            dlink = "ftp://data.argo.org.cn/pub/ARGO/BOA_Argo/NetCDF/"*f
            append!(tlist,[(dlink,settings["dir_analysis"]*f)])
        end
    end
    asyncmap(dl_BOA_file,tlist,ntasks=Threads.nthreads())
    return nothing
end

function dl_BOA_file(t_elem)
    println("  Downloading "*t_elem[1]*"...")
    retry   = true
    n_tries = 0
    while retry
        try
            n_tries += 1
            Downloads.download(t_elem[1],t_elem[2])
            retry = false
        catch
            if n_tries < 10
                println(t_elem[1]*" FAILED, retrying...")
                retry = true;
            else
                println(t_elem[1]*" FAILED, tried 10 times")
                retry = false
            end
        end
    end
    return nothing
end

function convert_BOA(settings)
    grid = process_grid_BOA(settings)
    save_grid_info(grid,settings)
    convert_monthly_files_BOA(grid,settings)
end

function process_grid_BOA(settings)
    # Load file and compute all grid variables
    printstyled(" Processing grid properties...\n",color=:blue)
    fh = Dataset(settings["dir_analysis"]*"BOA_Argo_"*string(settings["tvec"][1,1])*"_"*lpad(settings["tvec"][1,2],2,"0")*".nc","r")
    # Lat,Lon
    θ = fh["lat"][:]
    ϕ = fh["lon"][:] 

    # Height and pressure
    p = repeat(reshape(fh["pres"][:],(1,1,:)),length(ϕ),length(θ),1);
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
    slm = isfinite.(nomissing(fh["temp"][:,:,:,1],NaN32))
    @. V[~slm] = 0;
    @. A[~slm[:,:,1]] = 0;
    @. Δz[~slm] = 0;
    @. Δp[~slm] = 0;

    # Store in structure
    grid=grid_struct(ϕ,θ,z,p,Δz,Δp,A,V,slm);
    close(fh)
    return grid
end

function convert_monthly_files_BOA(grid,settings)
    # Convert EN4 files to the generic format
    printstyled(" Converting to generic format...\n",color=:blue)
    @sync @distributed for tstep in axes(settings["tvec"],1)
        yr = settings["tvec"][tstep,1]
        mnth = settings["tvec"][tstep,2]
        println("  Converting year "*string(yr)*" month "*string(mnth)*"...")

        # Read file
        fh = Dataset(settings["dir_analysis"]*"BOA_Argo_"*string(yr)*"_"*lpad(mnth,2,"0")*".nc","r")
        CT = nomissing(fh["temp"][:,:,:,1],NaN32);
        SA = nomissing(fh["salt"][:,:,:,1],NaN32);
        close(fh)
            
        # Convert to SA,CT,ρ
        @. SA = grid.slm * gsw_sa_from_sp(SA,grid.p,$reshape(grid.ϕ,(:,1,1)),$reshape(grid.θ,(1,:,1)));
        @. CT = grid.slm * gsw_ct_from_t(SA,CT,grid.p);  
        ρ = @. grid.slm * convert(Float32,gsw_rho(SA,CT,grid.p));
        data_monthly = data_struct(SA,CT,ρ)
        save_analysis_fmt(settings["dir_analysis_fmt"],settings["model_name"],yr,mnth,grid,data_monthly) # Save file
    end
    return nothing
end