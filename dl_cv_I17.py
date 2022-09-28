# Copyright [2022] JPL
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

# ----------------------------------------------------------------------
# Download and convert the I17 temperature and salinity
# grids, based on Ishii et al. (2017) into the standard format that
# can be read by the other scripts.
#
# Ishii, M., Fukuda, Y., Hirahara, S., Yasui, S., Suzuki, T.,
# & Sato, K. (2017). Accuracy of Global Upper Ocean Heat Content
# Estimation Expected from Present Observational Data Sets.
# SOLA, 13(0), 163–167. https://doi.org/10.2151/sola.2017-030
#
# This data is updated annually, links need to be updated manually
# ----------------------------------------------------------------------
from netCDF4 import Dataset
import pygrib
import gsw
import numpy as np
import os
import multiprocessing as mp
import urllib.request
mp.set_start_method('fork')

def main():
    global settings
    def_settings()
    if settings['download_files']: download_files()
    if settings['convert_files']:
        save_grid()
        convert_files()
    return

def def_settings():
    global settings
    settings = {}
    settings['dir_analysis'] = os.getenv('HOME') + '/Data/Steric/I17/analysis/'  # EN4 source directory
    settings['dir_analysis_fmt'] = os.getenv('HOME') + '/Data/Steric/I17/analysis_fmt/'  # EN4 target directory
    settings['depth'] = np.array([0,10,20,30,50,75,100,125,150,200,250,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1750,2000])
    settings['years'] = np.arange(2005, 2022, 1)
    settings['nprocs'] = 8 # Number of processes

    settings['download_files'] = True # Download the raw source data
    settings['convert_files'] = True  # Convert the raw source data to monthly NetCDF files

    return

def download_files():
    global settings
    mp.Pool(settings['nprocs']).map(download_indiv_year, settings['years'])
    return

def download_indiv_year(year):
    print('Year '+str(year))
    global settings
    url="https://climate.mri-jma.go.jp/pub/ocean/ts/v7.3.1/"
    url_sal = url + "sal/grib2/sal." + str(year) + ".grb2"
    url_temp = url + "temp/grib2/temp." + str(year) + ".grb2"
    fn_sal = settings["dir_analysis"] + "sal." +str(year) + ".grb2"
    fn_temp = settings["dir_analysis"] + "temp." + str(year) + ".grb2"
    print('   Year '+str(year)+' downloading salinity...')
    null = urllib.request.urlretrieve(url_sal,fn_sal)
    print('   Year '+str(year)+' downloading temperature...')
    null = urllib.request.urlretrieve(url_temp,fn_temp)
    return


def convert_files():
    global settings
    mp.Pool(settings['nprocs']).map(convert_indiv_year, settings['years'])
    return


def convert_indiv_year(year):
    print('Year '+str(year))
    global settings
    # Grid paramaters
    depth = np.array([0,10,20,30,50,75,100,125,150,200,250,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1750,2000])
    lon = np.arange(0.5,360.5,1)
    lat = np.arange(-89.5,90.5,1)
    fn_sal = settings["dir_analysis"] + "sal." +str(year) + ".grb2"
    fn_temp = settings["dir_analysis"] + "temp." + str(year) + ".grb2"
    sal_data  = np.zeros([12, 26, 180, 360],dtype=np.float32)
    temp_data = np.zeros([12, 26, 180, 360],dtype=np.float32)
    print('   Year '+str(year)+' reading salinity...')
    grbs_sal = pygrib.open(fn_sal)
    for grb in grbs_sal:
        if grb.parameterNumber == 193:
            if grb.level < 2001:
                zlev = np.where(grb.level == depth)[0][0]
                tlev = grb.month - 1
                d_local = grb.values._get_data()
                d_local[d_local > 1e10] = np.nan
                sal_data[tlev, zlev, :, :] = d_local

    print('   Year '+str(year)+' reading temperature...')
    grbs_temp = pygrib.open(fn_temp)
    for grb in grbs_temp:
        if grb.parameterNumber == 192:
            if grb.level < 2001:
                zlev = np.where(grb.level == depth)[0][0]
                tlev = grb.month - 1
                d_local = grb.values._get_data()
                d_local[d_local > 1e10] = np.nan
                temp_data[tlev, zlev, :, :] = d_local

    # Convert to SA,CT
    print('   Year '+str(year)+' converting to SA/CT...')
    p = gsw.p_from_z(-depth[:,np.newaxis,np.newaxis],lat[np.newaxis,:,np.newaxis])
    rho_data = np.zeros(temp_data.shape)
    slm = temp_data[0,...] < 1000
    # Mask out Caspian Sea
    slm[:,125:140, 44:55] = False

    for mnth in np.arange(12):
        sal_data[mnth,...] = gsw.SA_from_SP(sal_data[mnth,...],p,lon[np.newaxis,np.newaxis,:],lat[np.newaxis,:,np.newaxis])
        temp_data[mnth,...] = gsw.CT_from_t(sal_data[mnth,...], temp_data[mnth,...],p)
        rho_data[mnth, ...] = gsw.rho(sal_data[mnth,...], temp_data[mnth,...],p)

        # sea-land mask to zero
        sal_data[mnth, ~slm] = 0
        temp_data[mnth, ~slm] = 0
        rho_data[mnth, ~slm] = 0

    # Save data in good format
    print('   Year '+str(year)+' saving data...')
    for mnth in np.arange(12):
        fn = settings['dir_analysis_fmt'] + 'I17_'+str(year)+"_"+str(mnth+1).zfill(2)+'.nc'
        fh = Dataset(fn,'w')
        fh.createDimension('lon', len(lon))
        fh.createDimension('lat', len(lat))
        fh.createDimension('lvl', len(depth))
        fh.createVariable('lon', 'f4', ('lon',), zlib=True, complevel=4)[:] = lon
        fh.createVariable('lat', 'f4', ('lat',), zlib=True, complevel=4)[:] = lat
        fh.createVariable('lvl', 'i4', ('lvl',), zlib=True, complevel=4)[:] = np.arange(1,len(depth)+1)
        fh.createVariable('SA', 'f4', ('lvl', 'lat', 'lon',), zlib=True, complevel=4)[:] = sal_data[mnth,...]
        fh.createVariable('CT', 'f4', ('lvl', 'lat', 'lon',), zlib=True, complevel=4)[:] = temp_data[mnth,...]
        fh.createVariable('rho', 'f4', ('lvl', 'lat', 'lon',), zlib=True, complevel=4)[:] = rho_data[mnth,...]
        fh.close()
    return


def save_grid():
    global settings
    depth = np.array([0,10,20,30,50,75,100,125,150,200,250,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1750,2000])
    lon = np.arange(0.5,360.5,1)
    lat = np.arange(-89.5,90.5,1)
    fn_temp = settings["dir_analysis"] + "temp." + str(settings['years'][0]) + ".grb2"

    grbs_temp = pygrib.open(fn_temp)
    temp_data = np.zeros([12, 26, 180, 360],dtype=np.float32)
    for grb in grbs_temp:
        if grb.parameterNumber == 192:
            if grb.level < 2001:
                zlev = np.where(grb.level == depth)[0][0]
                tlev = grb.month - 1
                d_local = grb.values._get_data()
                d_local[d_local > 1e10] = np.nan
                temp_data[tlev, zlev, :, :] = d_local

    slm = temp_data[0,...] < 1000
    slm[:,125:140, 44:55] = False # Caspian Sea

    depth_bnds = np.zeros([len(depth),2])
    depth_bnds[0,1] = 5
    for i in range(1,len(depth)-1):
        depth_bnds[i,0] = depth_bnds[i-1,1]
        depth_bnds[i,1] = depth[i]+(depth[i+1] - depth[i])/2
    depth_bnds[-1, 0] = depth_bnds[-2,1]
    depth_bnds[-1, 1] = 2000

    delta_z = np.diff(depth_bnds,axis=1)
    p_top = gsw.p_from_z(-depth_bnds[:,0,np.newaxis,np.newaxis],lat[np.newaxis,:,np.newaxis])
    p_bottom = gsw.p_from_z(-depth_bnds[:,1,np.newaxis,np.newaxis],lat[np.newaxis,:,np.newaxis])
    p = gsw.p_from_z(-depth[:,np.newaxis,np.newaxis],lat[np.newaxis,:,np.newaxis])
    delta_p = p_bottom - p_top

    delta_p = np.repeat(delta_p,len(lon),axis=2)
    delta_z = np.repeat(delta_z[:,:,np.newaxis],len(lon),axis=2)
    delta_z = np.repeat(delta_z,len(lat),axis=1)

    z = np.repeat(depth[:,np.newaxis,np.newaxis],len(lon),axis=2)
    z = np.repeat(z,len(lat),axis=1)
    p = np.repeat(p,len(lon),axis=2)

    area = grid_area(lat,lon)
    volume = area[np.newaxis,:,:] * delta_z

    area[~slm[0,:,:]] = 0
    volume[~slm] = 0
    delta_p[~slm] = 0
    delta_z[~slm] = 0

    fn = settings["dir_analysis_fmt"] + "I17_grid.nc"
    fh = Dataset(fn,"w")
    fh.createDimension('lon', len(lon))
    fh.createDimension('lat', len(lat))
    fh.createDimension('lvl', len(depth))

    fh.createVariable('lon', 'f4', ('lon',), zlib=True, complevel=4)[:] = lon
    fh.createVariable('lat', 'f4', ('lat',), zlib=True, complevel=4)[:] = lat
    fh.createVariable('lvl', 'i4', ('lvl',), zlib=True, complevel=4)[:] = np.arange(1, len(depth) + 1)

    fh.createVariable('z', 'f4', ('lvl','lat','lon',), zlib=True, complevel=4)[:] = z
    fh.createVariable('p', 'f4', ('lvl','lat','lon',), zlib=True, complevel=4)[:] = p

    fh.createVariable('Dz', 'f4', ('lvl','lat','lon',), zlib=True, complevel=4)[:] = delta_z
    fh.createVariable('Dp', 'f4', ('lvl','lat','lon',), zlib=True, complevel=4)[:] = delta_p

    fh.createVariable('A', 'f4', ('lat','lon',), zlib=True, complevel=4)[:] = area
    fh.createVariable('V', 'f4', ('lvl','lat','lon',), zlib=True, complevel=4)[:] = volume
    fh.createVariable('slm', 'i2', ('lvl','lat','lon',), zlib=True, complevel=4)[:] = slm
    fh.close()
    return

def grid_area(lat,lon):
    grid = np.mean(np.diff(lat))
    lon,lat = np.meshgrid(lon,lat)
    radius = 6371000
    area =  (np.deg2rad(lon+grid/2)-np.deg2rad(lon-grid/2)) * (np.sin(np.deg2rad(lat+grid/2)) - np.sin(np.deg2rad(lat-grid/2)))*radius**2
    return(area)

if __name__ == "__main__":
    main()