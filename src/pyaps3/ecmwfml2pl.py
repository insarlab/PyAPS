# **************************** LICENSE START ***********************************
#
# Copyright 2022 ECMWF. This software is distributed under the terms
# of the Apache License version 2.0. In applying this license, ECMWF does not
# waive the privileges and immunities granted to it by virtue of its status as
# an Intergovernmental Organization or submit itself to any jurisdiction.
#
# ***************************** LICENSE END ************************************


import xarray as xr
import numpy as np
import sys
import os

# The list of A and B coefficients defining the model levels
PV_COEFF_A = np.array([
    0.0, 2.000365, 3.102241, 4.666084, 6.827977, 9.746966, 13.605424, 18.608931, 24.985718, 32.98571,
    42.879242, 54.955463, 69.520576, 86.895882, 107.415741, 131.425507, 159.279404, 191.338562, 227.968948,
    269.539581, 316.420746, 368.982361, 427.592499, 492.616028, 564.413452, 643.339905, 729.744141, 823.967834,
    926.34491, 1037.201172, 1156.853638, 1285.610352, 1423.770142, 1571.622925, 1729.448975, 1897.519287,
    2076.095947, 2265.431641, 2465.770508, 2677.348145, 2900.391357, 3135.119385, 3381.743652, 3640.468262,
    3911.490479, 4194.930664, 4490.817383, 4799.149414, 5119.89502, 5452.990723, 5798.344727, 6156.074219,
    6526.946777, 6911.870605, 7311.869141, 7727.412109, 8159.354004, 8608.525391, 9076.400391, 9562.682617,
    10065.97852, 10584.63184, 11116.66211, 11660.06738, 12211.54785, 12766.87305, 13324.66895, 13881.33106,
    14432.13965, 14975.61523, 15508.25684, 16026.11523, 16527.32227, 17008.78906, 17467.61328, 17901.62109,
    18308.43359, 18685.71875, 19031.28906, 19343.51172, 19620.04297, 19859.39063, 20059.93164, 20219.66406,
    20337.86328, 20412.30859, 20442.07813, 20425.71875, 20361.81641, 20249.51172, 20087.08594, 19874.02539,
    19608.57227, 19290.22656, 18917.46094, 18489.70703, 18006.92578, 17471.83984, 16888.6875, 16262.04688,
    15596.69531, 14898.45313, 14173.32422, 13427.76953, 12668.25781, 11901.33984, 11133.30469, 10370.17578,
    9617.515625, 8880.453125, 8163.375, 7470.34375, 6804.421875, 6168.53125, 5564.382813, 4993.796875, 4457.375,
    3955.960938, 3489.234375, 3057.265625, 2659.140625, 2294.242188, 1961.5, 1659.476563, 1387.546875, 1143.25,
    926.507813, 734.992188, 568.0625, 424.414063, 302.476563, 202.484375, 122.101563, 62.78125, 22.835938,
    3.757813, 0.0, 0.0])

PV_COEFF_B = np.array([
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7e-06, 2.4e-05, 5.9e-05,
    0.000112, 0.000199, 0.00034, 0.000562, 0.00089, 0.001353, 0.001992, 0.002857, 0.003971, 0.005378, 0.007133,
    0.009261, 0.011806, 0.014816, 0.018318, 0.022355, 0.026964, 0.032176, 0.038026, 0.044548, 0.051773, 0.059728,
    0.068448, 0.077958, 0.088286, 0.099462, 0.111505, 0.124448, 0.138313, 0.153125, 0.16891, 0.185689, 0.203491,
    0.222333, 0.242244, 0.263242, 0.285354, 0.308598, 0.332939, 0.358254, 0.384363, 0.411125, 0.438391, 0.466003,
    0.4938, 0.521619, 0.549301, 0.576692, 0.603648, 0.630036, 0.655736, 0.680643, 0.704669, 0.727739, 0.749797,
    0.770798, 0.790717, 0.809536, 0.827256, 0.843881, 0.859432, 0.873929, 0.887408, 0.8999, 0.911448, 0.922096,
    0.931881, 0.94086, 0.949064, 0.95655, 0.963352, 0.969513, 0.975078, 0.980072, 0.984542, 0.9885, 0.991984, 0.995003,
    0.99763, 1.0])

def get_input_variable_list(fin):
    ds = xr.open_dataset(fin)
    var_list = list(ds.keys())
    ds = None
    var_list_unique = list(set(var_list))
    if 'lnsp' not in var_list_unique:
      print("Error - lnsp variable missing from input file -exiting")
      sys.exit()
    if len(var_list_unique) < 2:
      print("Error - Data variable missing from input file -exiting")
      sys.exit()
    return var_list_unique

def check_requested_levels(plevs):
    check_lev = True
    if len(plevs) > 1:
        error_msg = "Error - only specify 1 input pressure level to interpolate to"
    else:
        for lev in plevs:
           if lev < 0 or lev > 110000 :
              check_lev = False
              error_msg = "Error - negative values and large positive values for pressure are not allowed -exiting"
    if check_lev == False:
        print(error_msg)
        sys.exit()
    return check_lev

def vertical_interpolate(vcoord_data, interp_var, interp_level):
    """A function to interpolate sounding data from each station to
    every millibar. Assumes a log-linear relationship.

    Input
    -----
    vcoord_data : A 1D array of vertical level values (e.g. from ERA5 pressure at model levels at a point)
    interp_var : A 1D array of the variable to be interpolated to the  pressure level
    interp_level : A 1D array containing the vertical level to interpolate to

    Return
    ------
    interp_data : A 1D array that contains the interpolated variable on the interp_level
    """

    inp_vals = []
    for l in interp_level:
      # Set NaN for the pressure out of the range
      if l < np.min(vcoord_data):
        ip = np.NAN
      elif l > np.max(vcoord_data):
        ip = np.NAN
      else:
        # Make vertical coordinate data and grid level log variables
        lnp = np.log(vcoord_data)
        lnp_interval = np.log(l)
        # Use numpy to interpolate from observed levels to grid levels
        ip = np.interp(lnp_interval, lnp, interp_var)
      inp_vals.append(ip)
    return inp_vals

def get_ph_levs(sp, level):
    '''Return the presure at a given level and the next'''
    ph_lev = PV_COEFF_A[level - 1] + (PV_COEFF_B[level - 1] * sp)
    ph_levplusone = PV_COEFF_A[level] + (PV_COEFF_B[level] * sp)
    return ph_lev, ph_levplusone

def compute_z_level(lev, sp, t_level, q_level, z_h):
    '''Compute z at half & full level for the given level, based on t/q/sp'''
    # compute moist temperature
    t = (t_level * (1. + 0.609133 * q_level)).copy()

    # compute the pressures (on half-levels)
    ph_lev, ph_levplusone = get_ph_levs(sp, lev)

    if lev == 1:
        dlog_p = np.log(ph_levplusone / 0.1)
        alpha = np.log(2)
    else:
        dlog_p = np.log(ph_levplusone / ph_lev)
        alpha = 1. - ((ph_lev / (ph_levplusone - ph_lev)) * dlog_p)

    # R_D = 287.06
    t = t * 287.06

    # z_f is the geopotential of this full level
    # integrate from previous (lower) half-level z_h to the
    # full level
    z_f = z_h + (t * alpha)

    # z_h is the geopotential of 'half-levels'
    # integrate z_h to next half level
    z_h = z_h + (t * dlog_p)

    return z_h, z_f

def calculate_pressure_on_model_levels(sp):
    # Get the a and b coefficients from the pv array to calculate the model level pressure
    p_half=[]
    for i in range(len(PV_COEFF_A)):
        p_half.append(PV_COEFF_A[i] + PV_COEFF_B[i] * sp)
    p_ml=[]
    for hybrid in range(len(p_half) - 1):
        p_ml.append((p_half[hybrid + 1] + p_half[hybrid]) / 2.0)
    return np.stack(p_ml)


def calculate_interpolated_pressure_field(data_var_on_ml, data_p_on_ml,plevs):
    nlevs = len(data_var_on_ml.level)
    p_array = np.stack(data_p_on_ml, axis=2).flatten()

    # Flatten the data array to enable faster processing
    var_array = np.stack(data_var_on_ml, axis=2).flatten()
    no_grid_points =  int(len(var_array)/nlevs)
    nlats_values = data_var_on_ml.coords['latitude']
    nlons_values = data_var_on_ml.coords['longitude']
    nlats = len(nlats_values)
    nlons = len(nlons_values)

    # Iterate over the data, selecting one vertical profile at a time
    count = 0
    profile_count = 0
    interpolated_values=[]
    for point in range(no_grid_points):
        offset =  count*nlevs
        var_profile = var_array[offset:offset+nlevs]
        p_profile = p_array[offset:offset+nlevs]
        interpolated_values.append(vertical_interpolate(p_profile, var_profile, plevs))
        profile_count += len(p_profile)
        count = count + 1
    interpolated_field=np.asarray(interpolated_values).reshape(nlats,nlons,len(plevs))
    interpolated_field = np.moveaxis(interpolated_field, -1, 0)
    return interpolated_field

def check_data_cube(dc):
    checks = True
    for var_name in list(dc.keys()):
        if var_name in ['z', 't', 'q', 'lnsp']:
            lnsp_dims = ['time', 'level', 'latitude','longitude']
            if all(value in lnsp_dims for value in dc[var_name].dims):
                continue
            else:
                print("Not all required lnsp dimensions found -exiting ",
                      dc[var_name].dims)
                checks = False
    return checks

def ECMWFml2pl(input_fname, output_fname, plevels):
    '''
    Convert the ECMWF model level netcdf to pressure level netcdf

    Args:
        * input_fname   : input ECMWF file at model level in netcdf format (str)
        * output_fname  : output pressure level in netcdf format (str)
        * plevels       : pressure levels in Pa (list)
    '''

    check_requested_levels(plevels)

    if not os.path.isfile(input_fname):
        print("Input file does not exist - exiting")
        sys.exit()

    variable_list = get_input_variable_list(input_fname)

    # Create a data object to hold the input and output data
    data_cube = xr.open_dataset(input_fname)

    if not check_data_cube(data_cube):
        sys.exit()
   # Get the ln surface pressure
    data_cube['pml']=data_cube['lnsp'].copy()
    data_cube['gph'] = xr.DataArray(
       data=np.zeros(data_cube['pml'].shape), dims=data_cube['pml'].dims,
       coords=data_cube['pml'].coords)

    for time_step in range(len(data_cube['t'].time)):
        data_slice_lnsp=data_cube['lnsp'][time_step,0,:,:]
        z_h=data_cube['z'][time_step,0,:,:]
        sp = np.exp(data_slice_lnsp)
        #  Get the pressure field on model levels for each timestep
        data_cube['pml'][time_step,:,:,:] = calculate_pressure_on_model_levels(sp)
        for lev in range(len(data_cube['t'].level), 0, -1):
            t_level = data_cube['t'][time_step,lev - 1,:,:]
            q_level = data_cube['q'][time_step,lev - 1,:,:]
            z_h, data_cube['gph'][time_step, lev - 1,:,:] = \
                compute_z_level(lev, sp, t_level, q_level, z_h)

    data_cube['pml'].attrs = {'units' : 'Pa','long_name':'pressure',
                              'standard_name':'air_pressure','positive':'down'}
    data_cube['gph'].attrs = {'units' : 'm**2 s**-2','long_name':'Geopotential',
                              'standard_name':'geopotential'}

    all_interpolated_var_fields_list = []
    for var in variable_list + ['gph']:
        if var in ['lnsp', 'pml', 'z']:
            continue
        interpolated_var_field = data_cube[var].copy()
        interpolated_var_field = interpolated_var_field[:,0:len(plevels),:,:]
        interpolated_var_field = interpolated_var_field.rename({'level':'pressure'})
        interpolated_var_field['pressure'] = plevels
        for time_step in range(len(data_cube[var].time)):
            var_on_ml = data_cube[var][time_step,:,:,:]
            p_on_ml = data_cube['pml'][time_step,:,:,:]
            interpolated_var_field[time_step,:,:,:] = \
                calculate_interpolated_pressure_field(var_on_ml,p_on_ml,plevels)
        all_interpolated_var_fields_list.append(interpolated_var_field)
    all_interpolated_var_fields = xr.merge(all_interpolated_var_fields_list)
    all_interpolated_var_fields = all_interpolated_var_fields.rename({'gph':'z'})
    all_interpolated_var_fields['pressure'].attrs = {'units' : 'Pa','long_name':'pressure',
                                                     'standard_name':'air_pressure','positive':'down'}
    all_interpolated_var_fields.to_netcdf(output_fname,format='NETCDF', engine='scipy',
                                          encoding={'z': {'zlib': True, 'complevel': 4}})
    all_interpolated_var_fields.close()
