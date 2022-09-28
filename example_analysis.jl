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

# Simple example to plot and analyze the output data
using GMT
using NCDatasets
using LinearAlgebra
using Statistics
using Dates

# Open the file and read the data
fh = Dataset(homedir()*"/Data/Steric/SIO/postprocessed/SIO_steric_ohc_2004_2021.nc","r")
t = nomissing(fh["time"][:])
ζ = nomissing(fh["ohc_ts"][:])
ϕ = nomissing(fh["lon"][:])
θ = nomissing(fh["lat"][:])
η = nomissing(fh["totalsteric_2d"][:])
close(fh)

# Plot global-mean ocean heat content
plot(t, ζ/1e23, lw=1.5, lc="#4c78a8", ylabel="Global ocean heat content (ZJ)", show=true, limits=(t[1],t[end],-1.2,1.4),figname="ohc_timeseries.png",figsize=(5,4),conf=(FONT_ANNOT=6,FONT_LABEL=6,PS_LINE_JOIN=round))

# Compute the linear trend in steric sea level.
tval = @. year(t) + month(t)/12 - 1/24
design_matrix = ones(length(tval),6)
@. design_matrix[:,2] = tval - $mean(tval)
@. design_matrix[:,3] = sin(2π*tval)
@. design_matrix[:,4] = cos(2π*tval)
@. design_matrix[:,5] = sin(4π*tval)
@. design_matrix[:,6] = cos(4π*tval)
qrd = qr(design_matrix)
rinv = inv(qrd.R) * qrd.Q'
steric_trend = zeros(length(ϕ),length(θ))
for i in CartesianIndices(steric_trend)
    steric_trend[i] = (rinv * η[i,:])[2]
end

# Plot map with GMT
cpt = makecpt(color=:vik, range=(-5,5,0.25))
grid = mat2grid(steric_trend'.*1000,x=convert.(Float64,ϕ),y=convert.(Float64,θ))
grdimage(grid,projection=:robinson,color=cpt,region=(-180, 180, -90, 90),figsize=(6.5,4),conf=(FONT_ANNOT_PRIMARY=6,FONT_LABEL=6))
coast!(land=:darkgrey)
colorbar!(pos=(anchor=:BC,length=(4,0.2), horizontal=true, offset=(0,0.5),triangles=true), color=cpt, frame=(annot=2, xlabel="Steric sea level (mm yr@+-1@+)"),  show=true, nolines=true,figname="trend_map.png",conf=(FONT_ANNOT_PRIMARY=6,FONT_LABEL=6))