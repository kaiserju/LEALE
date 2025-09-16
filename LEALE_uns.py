'''
Created by Júlia Kaiser & Fabio Viola - Jul/2025
CMCC Foundation - Euro-Mediterranean Center on Climate Change
GOCO - Global Coastal Ocean Division
'''
import os
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyresample import geometry, kd_tree
from geopy.distance import geodesic
from datetime import datetime 

def haversine(lat1, lon1, lat2, lon2):
    R = 6371000  # Earth radius in meters
    phi1, phi2 = np.radians(lat1), np.radians(lat2)
    dphi = phi2 - phi1
    dlambda = np.radians(lon2 - lon1)

    a = np.sin(dphi/2)**2 + np.cos(phi1)*np.cos(phi2)*np.sin(dlambda/2)**2
    return 2 * R * np.arcsin(np.sqrt(a))

#############################
## FOLDERS + EXPERIMENT ID ##
#############################
pathnamein = input('Pathname of unstructured input .nc files: ')
pathnameout = input('Pathname to save LEALE results: ')
pathnameETICOfile = input('Complete filename of thalweg (from ETICO tool): ')

files = os.listdir(pathnamein)
nosfiles = []
for file in files:
    if file.endswith('.nc'):
        nosfiles.append(file)
nosfiles.sort()


########################
## OPPENING VARIABLES ##
########################
vec_time = [];
vec_sal = np.array([]);
vec_lat = []; vec_lon = []; vec_distance = []

## Choosing the salinity threshold to consider as SWI
SWI_threshold = float(input('Enter with salinity threshold to be considered: '))
#SWI_threshold = 1.0


############################
## FINDING NEAREST POINTS ##
##### UNSTRUCTURED GRID ####
############################

# Loading thalweg from ETICO tool (the points should be ordered from the river's mouth to the head!)
rightorder=True #CHANGE HERE WHEN NEEDED!
thalweg = pd.read_csv(pathnameETICOfile)
lat_etico_extr = np.array(thalweg['Latitude'])
lon_etico_extr = np.array(thalweg['Longitude'])

if rightorder==True:
    lat_etico = lat_etico_extr
    lon_etico = lon_etico_extr
else:
    lat_etico = lat_etico_extr[::-1]
    lon_etico = lon_etico_extr[::-1]

# River mouth coordinates
river_mouth_coords_etico = (lat_etico[0], lon_etico[0])

# Loading model result (unstructured grid)
uns = nc.Dataset(pathnamein + nosfiles[0]) #CHANGE HERE iF NEEDED!
lat_uns = uns.variables['latitude'][:].data
lon_uns = uns.variables['longitude'][:].data
levels = uns.variables['level']
depth = uns.variables['total_depth']
del uns

# Define the swath (unstructured grid)
swath_def = geometry.SwathDefinition(lons=lon_uns, lats=lat_uns)

# Get the nearest neighbors
target_geo = geometry.SwathDefinition(lons=lon_etico, lats=lat_etico)
_, _, indices, distances = kd_tree.get_neighbour_info(swath_def, target_geo, radius_of_influence=30, neighbours=1)
indices = np.array(indices, dtype=int)

vec_sal = np.zeros((0,len(indices),0)) #[time,point,sal]
depth_indices = depth[indices][:].data
lastlevelixd_indices = []
for i in range(0,len(depth_indices)):
    idx = np.where(levels[:].data <= depth_indices[i])[0][-1] #TAKING THE LAST ONE (BOTTOM)
    lastlevelixd_indices.append(idx)


## Getting the distance between each point and the point I'm considering as the river mouth
distance_from_rivermouth = []; index_thalweg = []
for d in range(0,len(lat_etico)):
    if d == 0:
        current_tmpdist = haversine(river_mouth_coords_etico[0], river_mouth_coords_etico[1], lat_uns[indices[d]], lon_uns[indices[d]])
        previous_tmpdist = haversine(river_mouth_coords_etico[0], river_mouth_coords_etico[1], lat_uns[indices[d]], lon_uns[indices[d]])
    else:
        previous_tmpdist = distance_from_rivermouth[d-1]
        current_tmpdist = haversine(lat_uns[indices[d-1]], lon_uns[indices[d-1]], lat_uns[indices[d]], lon_uns[indices[d]])
    tmpdist = current_tmpdist + previous_tmpdist
    distance_from_rivermouth.append(tmpdist)
    index_thalweg.append(d)


###########################################################################
############### COMPUTING SWI length FOR ALL SHYFEM OUTPUT ################
### METHOD CONSIDERING BOTH (THE BOTTOM LAYER + HIGHEST SALINITY LAYER) ###
###########################################################################
maxsal_belowthreshold = []; bottomsal_belowthreshold = []
results = []

for file in nosfiles:
    # Loading files
    uns = nc.Dataset(pathnamein + file)
    salinity = uns.variables['salinity'][:].data #[time,points,depth]
    time_data = uns.variables['time']
    levels = uns.variables['level'][:].data
    levels = np.round(levels[:].data,2)
    depth = uns.variables['total_depth'][:].data
    latitude = uns.variables['latitude'][:].data
    longitude = uns.variables['longitude'][:].data
    del uns
    
    ii = 0
    for p in indices: #loop through the thalweg points
        salinity_point = salinity[:, p, :]
        depth_point = np.round(depth[p], 2)
        last_active_layer = np.max(np.where(levels <= depth_point))
        
        if (np.max(salinity[:,p,last_active_layer+1]) != 0) & ((float(levels[last_active_layer+1]) - depth_point) < 0.5): #a better way to certify we are getting the last active layer
            last_active_layer = last_active_layer + 1

        max_salinity_layer_idx = np.argmax(salinity_point, axis=1)  #max layer index over time
        max_salinity_value = np.max(salinity_point, axis=1)  #max salinity over time
        salinity_last_active_layer = salinity_point[np.arange(salinity_point.shape[0]), last_active_layer]  #value at last active layer

        # Calculate whether salinity is below the threshold
        maxsal_belowthreshold = (max_salinity_value < SWI_threshold).astype(int)
        bottomsal_belowthreshold = (salinity_last_active_layer < SWI_threshold).astype(int) #0 = sal. > threshold; 1 = sal. < threshold
        
        
        for t in range(0, len(time_data)):
            time_datetime = nc.num2date(time_data[t], units=time_data.units, calendar=time_data.calendar)
            time_datetime = datetime.strptime(str(time_datetime), '%Y-%m-%d %H:%M:%S')
            
                # Append results for this point and time step
            results.append({
                'latitude': latitude[p],
                'longitude': longitude[p],
                'time': time_datetime,
                'last_active_layer_idx': last_active_layer,
                'max_salinity_layer_idx': max_salinity_layer_idx[t],
                'max_salinity_value': max_salinity_value[t],
                'salinity_last_active_layer': salinity_last_active_layer[t],
                'max_salinity_below_SWIthreshold': maxsal_belowthreshold[t],
                'bottom_salinity_below_SWIthreshold': bottomsal_belowthreshold[t],
                'distance_from_river_mouth': distance_from_rivermouth[ii],
                'index_thalweg': index_thalweg[ii],               
            })
            
        ii=ii+1
            
        try:
            print(int(np.where(indices == p)[0]),  ' --- ', len(indices))
        except:
            print(int(np.where(indices == p)[0][0]), ' --- ', len(indices))
df = pd.DataFrame(results)

######################################################
## CLEANING THE DATAFRAME AND CHOOSING INFO TO SAVE ##
######################################################
df_results = df.copy()
times = df_results['time'].unique()
    
df_bottom = pd.DataFrame()
df_maxsallayer = pd.DataFrame()
for tt in times:
    #taking all information from a single time step
    df_tmp = df_results.loc[df_results['time'] == tt]
    
    ## --- Bottom layer --- ##
    #taking only the points at this time step that are with salinity higher the threshold (1PSU)
    df_tmp1 = df_tmp.loc[df_tmp['bottom_salinity_below_SWIthreshold'] == 0]
    if np.any(df_tmp1):
        df_tmp1 = df_tmp1.loc[df_tmp1['distance_from_river_mouth'] == np.max(df_tmp1['distance_from_river_mouth'])][0:1]
    
    else:
        print('df_tmp1 is empty! No elements with sal. > 1 were found along the thalweg path.')
        df_tmp1 = pd.DataFrame(columns=['latitude', 'longitude', 'time', 'last_active_layer_idx',
                                        'max_salinity_layer_idx', 'max_salinity_value',
                                        'salinity_last_active_layer', 'max_salinity_below_SWIthreshold',
                                        'bottom_salinity_below_SWIthreshold', 'distance_from_river_mouth',
                                        'index_thalweg'])
        df_tmp1.loc[0] = np.nan
        df_tmp1['time'][0] = tt
        df_tmp1['distance_from_river_mouth'] = 0
        
    if np.any(df_bottom) == False:
        df_bottom = df_tmp1.copy()
    else:
        df_bottom = pd.concat([df_bottom, df_tmp1])       
        
     ## --- Highest salinity layer --- ##
    #taking only the points at this time step that are with salinity higher the threshold (1PSU)
    df_tmp2 = df_tmp.loc[df_tmp['max_salinity_below_SWIthreshold'] == 1]
    if np.any(df_tmp2):
        #df_tmp2['max_salinity_value'] = np.round(df_tmp2['max_salinity_value'], 2)
        #df_tmp2 = df_tmp2.loc[df_tmp2['salinity_last_active_layer'] > 0]
        #taking the main SWI front
        #df_tmp2 = df_tmp2.loc[df_tmp2['distance_from_river_mouth'] == np.min(df_tmp2['distance_from_river_mouth'])][0:1]
        df_tmp2 = df_tmp2.loc[df_tmp2['distance_from_river_mouth'] == np.max(df_tmp2['distance_from_river_mouth'])][0:1]
    
    else:
        print('df_tmp2 is empty! No elements with sal. > 1 were found along the thalweg path.')
        df_tmp2 = pd.DataFrame(columns=['latitude', 'longitude', 'time', 'last_active_layer_idx',
       'max_salinity_layer_idx', 'max_salinity_value',
       'salinity_last_active_layer', 'max_salinity_below_SWIthreshold',
       'bottom_salinity_below_SWIthreshold', 'distance_from_river_mouth',
       'index_thalweg', 'distance_from_river_mouth_FV'])
        df_tmp2.loc[0] = np.nan
        df_tmp2['time'][0] = tt
        df_tmp2['distance_from_river_mouth'] = 0
        
    if np.any(df_maxsallayer) == False:
        df_maxsallayer = df_tmp2.copy()
    else:
        df_maxsallayer = pd.concat([df_maxsallayer, df_tmp2])

#df_bottom = df_bottom[['latitude','longitude','time','last_active_layer_idx','max_salinity_layer_idx','max_salinity_value','salinity_last_active_layer','max_salinity_below_SWIthreshold','bottom_salinity_below_SWIthreshold','distance_from_river_mouth','distance_from_river_mouth_FV','index_thalweg']]
#df_maxsallayer = df_maxsallayer[['latitude','longitude','time','last_active_layer_idx','max_salinity_layer_idx','max_salinity_value','salinity_last_active_layer','max_salinity_below_SWIthreshold','bottom_salinity_below_SWIthreshold','distance_from_river_mouth','distance_from_river_mouth_FV','index_thalweg']]

# Saving file
df_bottom = df_bottom[['latitude','longitude','time','salinity_last_active_layer','bottom_salinity_below_SWIthreshold','distance_from_river_mouth']]
df_bottom = df_bottom.loc[df_bottom['bottom_salinity_below_SWIthreshold'] == 0]; del df_bottom['bottom_salinity_below_SWIthreshold']
df_bottom.columns = ['Latitude','Longitude','Time','Salinity [PSU]','SWI length [m]']; df_bottom = df_bottom.reset_index(); del df_bottom['index']
df_bottom.to_csv(pathnameout + 'SWIleght_from_LEALE.csv', index=False)
