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

function download_EN4(settings)
    printstyled(" Downloading EN4 data...\n",color=:blue)

    # Download list with files:
    settings["model_name"] == "EN4_l09" ? flist = "https://www.metoffice.gov.uk/hadobs/en4/EN.4.2.2.analyses.l09.download-list.txt" : nothing
    settings["model_name"] == "EN4_g10" ? flist = "https://www.metoffice.gov.uk/hadobs/en4/EN.4.2.2.analyses.g10.download-list.txt" : nothing
    settings["model_name"] == "EN4_c13" ? flist = "https://www.metoffice.gov.uk/hadobs/en4/EN.4.2.2.analyses.c13.download-list.txt" : nothing
    settings["model_name"] == "EN4_c14" ? flist = "https://www.metoffice.gov.uk/hadobs/en4/EN.4.2.2.analyses.c14.download-list.txt" : nothing
    Downloads.download(flist,settings["dir_analysis"]*"dlist.txt")
    dlist = readdlm(settings["dir_analysis"]*"dlist.txt")
    acclist = zeros(Bool,length(dlist))
    for i in eachindex(acclist)
        yr = parse(Int32,dlist[i][end-7:end-4])
        acclist[i] = (yr in settings["years"])
    end
    dlist = dlist[acclist]

    @sync @distributed for idx in eachindex(dlist)
        yr = parse(Int32,dlist[idx][end-7:end-4])
        println("   Downloading year "*string(yr)*"...")
        url = convert(String,dlist[idx])
        svf = settings["dir_analysis"]*settings["model_name"]*"_"*string(yr)*".zip"
        Downloads.download(url,svf)
        zip_handle = ZipFile.Reader(svf);
        for f in zip_handle.files
            write(settings["dir_analysis"]*f.name,read(f, String));
        end
        close(zip_handle)
        rm(svf)
    end
    printstyled(" Downloading EN4 data done\n",color=:blue)
    return nothing
end

function convert_EN4(settings)
    xbt = settings["model_name"][end-2:end]
    grid = process_grid_EN4(xbt,settings)
    save_grid_info(grid,settings)
    convert_monthly_files_EN4(xbt,grid,settings)
end

function process_grid_EN4(xbt,settings)
    # Load file and compute all grid variables
    printstyled(" Processing grid properties...\n",color=:blue)
    fh = Dataset(settings["dir_analysis"]*"EN.4.2.2.f.analysis."*xbt*"."*string(settings["years"][1])*"01.nc","r")
    θ = fh["lat"][:]
    ϕ = fh["lon"][:]

    # Height and pressure
    z = repeat(reshape(fh["depth"][1:30],(1,1,:)),length(ϕ),length(θ),1);
    p = @. convert(Float32,gsw_p_from_z(-z, $reshape(θ,(1,:,1))));
    
    # Height and pressure levels
    z_bnds = fh["depth_bnds"][:,1:30];
    Δz = repeat(reshape(diff(z_bnds,dims=1)[1,:],(1,1,:)),length(ϕ),length(θ),1);

    ptop = repeat((@. gsw_p_from_z($reshape(-z_bnds[1,:],(1,1,:)), $reshape(θ,(1,:,1)))), length(ϕ),1,1);
    pbot = repeat((@. gsw_p_from_z($reshape(-z_bnds[2,:],(1,1,:)), $reshape(θ,(1,:,1)))), length(ϕ),1,1);
    Δp = @. convert(Float32,pbot - ptop);

    # Area and volume
    A = grid_area(ϕ,θ);
    V = @. A * Δz;

    # Sea-land mask
    slm = ismissing.(fh["temperature"][:,:,1:30,1]) .== 0 ;
    @. V[~slm] = 0;
    @. A[~slm[:,:,1]] = 0;
    @. Δz[~slm] = 0;
    @. Δp[~slm] = 0;

    # Store in structure
    grid=grid_struct(ϕ,θ,z,p,Δz,Δp,A,V,slm);
    return grid
end

function convert_monthly_files_EN4(xbt,grid,settings)
    # Convert EN4 files to the generic format
    printstyled(" Converting to generic format...\n",color=:blue)
    @sync @distributed for tstep in axes(settings["tvec"],1)
        yr = settings["tvec"][tstep,1]
        mnth = settings["tvec"][tstep,2]
        println("  Converting year "*string(yr)*" month "*string(mnth)*"...")

        # Read file
        fh = Dataset(settings["dir_analysis"] * "EN.4.2.2.f.analysis."*xbt*"."*string(yr)*lpad(string(mnth),2,"0")*".nc","r")
        CT = nomissing(fh["temperature"][:,:,1:30,1],0).-273.15
        SA = nomissing(fh["salinity"][:,:,1:30,1],0)
        close(fh)

        # Convert to SA,CT,ρ
        @. SA = grid.slm * gsw_sa_from_sp(SA,grid.p,$reshape(grid.ϕ,(:,1,1)),$reshape(grid.θ,(1,:,1)));
        @. CT = grid.slm * gsw_ct_from_pt(SA,CT);  
        ρ = @. grid.slm * convert(Float32,gsw_rho(SA,CT,grid.p));
        data_monthly = data_struct(SA,CT,ρ)
        save_analysis_fmt(settings["dir_analysis_fmt"],settings["model_name"],yr,mnth,grid,data_monthly) # Save file
    end
    return nothing
end

