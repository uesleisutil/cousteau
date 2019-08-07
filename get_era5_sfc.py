import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type':'reanalysis',
        'variable':[
            '10m_u_component_of_wind','10m_v_component_of_wind','2m_dewpoint_temperature',
            '2m_temperature','convective_snowfall','convective_snowfall_rate_water_equivalent',
            'friction_velocity','ice_temperature_layer_1','ice_temperature_layer_2',
            'ice_temperature_layer_3','ice_temperature_layer_4','land_sea_mask',
            'large_scale_snowfall','large_scale_snowfall_rate_water_equivalent','mean_sea_level_pressure',
            'sea_ice_cover','sea_surface_temperature','skin_temperature',
            'snow_albedo','snow_density','snow_depth',
            'snow_evaporation','snowfall','snowmelt',
            'soil_temperature_level_1','soil_temperature_level_2','soil_temperature_level_3',
            'soil_temperature_level_4','soil_type','surface_pressure',
            'temperature_of_snow_layer','total_column_snow_water','total_precipitation',
            'volumetric_soil_water_layer_1','volumetric_soil_water_layer_2','volumetric_soil_water_layer_3',
            'volumetric_soil_water_layer_4'
        ],
        'year':'2015',
        'month':[
            '10','11'
        ],
        'day':[
            '01','02','03',
            '04','05','06',
            '07','08','09',
            '10','11','12',
            '13','14','15',
            '16','17','18',
            '19','20','21',
            '22','23','24',
            '25','26','27',
            '28','29','30',
            '31'
        ],
        'time':[
            '00:00','06:00','12:00',
            '18:00'
        ],
        'format':'grib'
    },
    'download.grib')
