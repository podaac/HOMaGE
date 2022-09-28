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

# ------------------------------------------------------------------------
# Download and convert the `Institute of Atmospheric Physics, 
# Chinese Academy of Sciences` temperature and salinity grids
# into the standard format that can be read by the other scripts
#
# Cheng, L., Trenberth, K. E., Gruber, N., Abraham, J. P., Fasullo, 
# J. T., Li, G., Mann, M. E., Zhao, X., & Zhu, J. (2020). Improved 
# Estimates of Changes in Upper Ocean Salinity and the Hydrological 
# Cycle. Journal of Climate, 33(23), 10357–10381. 
# https://doi.org/10.1175/JCLI-D-20-0366.1
# 
# Cheng, L., & Zhu, J. (2016). Benefits of CMIP5 Multimodel Ensemble 
# in Reconstructing Historical Ocean Subsurface Temperature Variations. 
# Journal of Climate, 29(15), 5393–5416.
# https://doi.org/10.1175/JCLI-D-15-0730.1
# ------------------------------------------------------------------------

function download_IAP(settings)
    printstyled(" Downloading IAP data...\n",color=:blue)
    @sync @distributed for t in axes(settings["tvec"],1)
        println("   Downloading year "*string(settings["tvec"][t,1])*" month "*string(settings["tvec"][t,2])*" ...")
        fn_t = settings["dir_analysis"]*"CZ16_1_2000m_Temp_year_"*string(settings["tvec"][t,1])*"_month_"*lpad(settings["tvec"][t,2],2,"0")*".nc"
        url_t = "http://www.ocean.iap.ac.cn/ftp/cheng/CZ16_v3_IAP_Temperature_gridded_1month_netcdf/Monthly/CZ16_1_2000m_Temp_year_"*string(settings["tvec"][t,1])*"_month_"*lpad(settings["tvec"][t,2],2,"0")*".nc"
        Downloads.download(url_t,fn_t)
        fn_s = settings["dir_analysis"]*"CZ16_1_2000m_salinity_year_"*string(settings["tvec"][t,1])*"_month_"*lpad(settings["tvec"][t,2],2,"0")*".nc"
        url_s = "http://www.ocean.iap.ac.cn/ftp/cheng/CZ16_v0_IAP_Salinity_gridded_1month_netcdf/Monthly/CZ16_1_2000m_salinity_year_"*string(settings["tvec"][t,1])*"_month_"*lpad(settings["tvec"][t,2],2,"0")*".nc"
        Downloads.download(url_s,fn_s)
    end
    return nothing
end

function convert_IAP(settings)
    grid = process_grid_IAP(settings)
    save_grid_info(grid,settings)
    convert_monthly_files_IAP(grid,settings)
end

function process_grid_IAP(settings)
    # Load file and compute all grid variables
    printstyled(" Processing grid properties...\n",color=:blue)
    fh_t = Dataset(settings["dir_analysis"]*"CZ16_1_2000m_Temp_year_"*string(settings["tvec"][1,1])*"_month_"*lpad(settings["tvec"][1,2],2,"0")*".nc","r")
    fh_s = Dataset(settings["dir_analysis"]*"CZ16_1_2000m_salinity_year_"*string(settings["tvec"][1,1])*"_month_"*lpad(settings["tvec"][1,2],2,"0")*".nc","r")
    θ = nomissing(fh_s["lat"][:])
    ϕ = nomissing(fh_s["lon"][:])
    z = fh_s["depth_std"][:]
    Δz = zeros(Float32,2,length(z))
    Δz[2,1] = (z[2] - z[1])/2
    for k in 2:length(z)-1
        Δz[1,k] = Δz[2,k-1]
        Δz[2,k] = z[k] +(z[k+1] - z[k])/2
    end
    Δz[1,end] = Δz[2,end-1]
    Δz[2,end] = 2000
    ptop = @. gsw_p_from_z(-$reshape(Δz[1,:],(1,1,:)), $reshape(θ,(1,:,1)));
    pbot = @. gsw_p_from_z(-$reshape(Δz[2,:],(1,1,:)), $reshape(θ,(1,:,1)));
    Δp = repeat(pbot - ptop,length(ϕ),1,1);
    Δz = repeat(reshape(diff(Δz,dims=1)[1,:],(1,1,:)),length(ϕ),length(θ),1);
    p = repeat((@. gsw_p_from_z(-$reshape(z,(1,1,:)), $reshape(θ,(1,:,1)))),length(ϕ),1,1);;
    A = grid_area(ϕ,θ);
    V = @. $reshape(A,(360,180,1)) * Δz;
    slm = isfinite.(flip_field(nomissing(fh_s["salinity"][:],NaN32) .+ nomissing(fh_t["temp"][:],NaN32)))
    z = repeat(reshape(z,(1,1,:)),length(ϕ),length(θ),1);
    @. V[~slm] = 0;
    @. Δz[~slm] = 0;
    @. Δp[~slm] = 0;
    @. A[~slm[:,:,1]] = 0;
    # Depth
    dpth = zeros(Float32,size(A))
    for i=1:length(ϕ),j=1:length(θ)
        if slm[i,j,1]
            dpth[i,j] = z[findall(slm[i,j,:] )[end]]
        end
    end    
    grid=grid_struct(ϕ,θ,z,p,Δz,Δp,A,V,slm);
    close(fh_t)
    close(fh_s)
    return grid
end

function convert_monthly_files_IAP(grid,settings)
    # Convert EN4 files to the generic format
    printstyled(" Converting to generic format...\n",color=:blue)
    @sync @distributed for t in axes(settings["tvec"],1)
        yr = settings["tvec"][t,1]
        mnth = settings["tvec"][t,2]
        println("  Converting year "*string(yr)*" month "*string(mnth)*"...")
        fh_t = Dataset(settings["dir_analysis"]*"CZ16_1_2000m_Temp_year_"*string(settings["tvec"][t,1])*"_month_"*lpad(settings["tvec"][t,2],2,"0")*".nc","r")
        fh_s = Dataset(settings["dir_analysis"]*"CZ16_1_2000m_salinity_year_"*string(settings["tvec"][t,1])*"_month_"*lpad(settings["tvec"][t,2],2,"0")*".nc","r")

        # Read file
        SA = flip_field(nomissing(fh_s["salinity"][:],NaN32));
        CT = flip_field(nomissing(fh_t["temp"][:],NaN32));
        @. SA[~grid.slm] = 0
        @. CT[~grid.slm] = 0
    
        # Convert to SA,CT,ρ
        @. SA = grid.slm * gsw_sa_from_sp(SA,grid.p,$reshape(grid.ϕ,(:,1,1)),$reshape(grid.θ,(1,:,1)));
        @. CT = grid.slm * gsw_ct_from_t(SA,CT,grid.p);  
        ρ = @. grid.slm * convert(Float32,gsw_rho(SA,CT,grid.p));
        data_monthly = data_struct(SA,CT,ρ)
        save_data_file(mnth,yr,settings,grid,data_monthly) # Save file
    end
    return nothing
end

@everywhere flip_field = (field) -> permutedims(field, [2, 3, 1])