# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import KratosMultiphysics
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
from KratosMultiphysics.MeshMovingApplication.mesh_moving_analysis import MeshMovingAnalysis
import KratosMultiphysics.StatisticsApplication

import sys
import time
from os import makedirs as os_makedirs
import math


class FluidDynamicsAnalysisWithFlush(FluidDynamicsAnalysis):

    def __init__(self, model, project_parameters, flush_frequency=10.0):
        super(FluidDynamicsAnalysisWithFlush, self).__init__(
            model, project_parameters)
        self.flush_frequency = flush_frequency
        self.last_flush = time.time()

    def FinalizeSolutionStep(self):
        super(FluidDynamicsAnalysisWithFlush, self).FinalizeSolutionStep()

        if self.parallel_type == "OpenMP":
            now = time.time()
            if now - self.last_flush > self.flush_frequency:
                sys.stdout.flush()
                self.last_flush = now


def ccw_rotate_point_around_z(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy


if __name__ == "__main__":
    from os import path

    ###############################################################
    #
    # Input data
    # initialize the static angle case setup

    # Permited cases
    # here using CCW (counterclockwise) RH (right-handed) cartesian coordinate system definition
    available_sa_cases = [0, 1, 2, 4, 6, 8, 10]
    available_u_cases = [7.63, 15.26]

    # Time parameters
    rampup_time = 0.82
    # total_time -> it is calculated later as 10*rampup_time

    # for the 2d test
    delta_t = 0.01
    # for the 3d test
    #delta_t = 0.003
    # for the 3d fine mesh
    #delta_t = 0.00015

    restart_time_label = '2.47'

    ###############################################################
    #
    # Check that inut arguments were given and are among the available cases

    try:
        static_case = float(sys.argv[1])
        U = float(sys.argv[2])
    except:
        msg = 'Static case needs to be provided as a float argument after ' + \
            path.basename(__file__) + '\n'
        msg += 'this can be one of: ' + \
            ', '.join([str(val) for val in available_sa_cases]) + '\n'
        msg += 'Velocity needs to be provided as a float argument after ' + \
            path.basename(__file__) + '\n'
        raise Exception(msg)

    if static_case not in available_sa_cases:
        msg = 'Static case ' + str(static_case) + ' not valid\n'
        msg += 'this has to be one of: ' + \
            ', '.join([str(val) for val in available_sa_cases])
        raise Exception(msg)

    if U not in available_u_cases:
        msg = 'U ' + str(U) + ' not valid\n'
        msg += 'this has to be one of: ' + \
            ', '.join([str(val) for val in available_u_cases])
        raise Exception(msg)

    if U == 7.63:
        rampup_time *= 2
        delta_t *= 2

    ###############################################################
    #
    # initialize the static angle case setup

    total_time = 10*rampup_time
    cfd_time = 8*rampup_time

    total_time_interval = [0, 10*rampup_time]
    ale_interval = [0, 2*rampup_time]
    cfd_interval = [2*rampup_time, 10*rampup_time]
    rampup_interval = [2*rampup_time, 3*rampup_time]
    constant_condition_interval = [2*rampup_time, 10*rampup_time]
    stabilisation_interval = [3*rampup_time, 4*rampup_time]
    converged_flow_interval = [4*rampup_time, 10*rampup_time]

    if static_case < 0:
        folder_ident = path.join('sa', 'angle_m' + str(abs(static_case)))
    elif static_case > 0:
        folder_ident = path.join('sa', 'angle_p' + str(static_case))
    else:
        folder_ident = path.join('sa', 'angle_' + str(static_case))

    msg = 'Configuring case setup:\n'
    msg += '    - Angle of attack: ' + str(static_case) + '°\n'
    msg += '    - ALE time interval: ' + str(ale_interval) + ' s\n'
    msg += '    - CFD time interval: ' + str(cfd_interval) + ' s\n'
    msg += '    - Rampup time interval: ' + str(rampup_interval) + ' s\n'
    msg += '    - Constant conditions time interval: ' + str(constant_condition_interval) + ' s\n'
    msg += '    - Stabilisation time interval: ' + str(stabilisation_interval) + ' s\n'
    msg += '    - Converged flow time interval: ' + str(converged_flow_interval) + ' s\n'

    KratosMultiphysics.Logger.PrintInfo('MainStaticAngles', msg)


    ###############################################################
    ###############################################################
    #
    # read and adapt mesh moving - ALE - project parameters
    with open("ProjectParametersOpenMP_SA_ALE.json", 'r') as parameter_file:
    #with open("ProjectParametersOpenMP_SA_ALERestart.json", 'r') as parameter_file:
        parameters_ale = KratosMultiphysics.Parameters(parameter_file.read())

    print('\n---------- ADAPTING ALE PARAMETERS ----------\n')

    # Changing time parameters
    old_end_time = parameters_ale['problem_data']['end_time'].GetDouble()
    old_delta_t = parameters_ale['solver_settings']['time_stepping']['time_step'].GetDouble()

    parameters_ale['problem_data']['end_time'].SetDouble(ale_interval[1])
    parameters_ale['solver_settings']['time_stepping']['time_step'].SetDouble(ale_interval[1])

    msg = 'Time parameters adapted:\n'
    msg += '    - ' + str(old_end_time) + ' was replaced by ' + str(ale_interval[1]) + ' as the end time.\n'
    msg += '    - ' + str(old_delta_t) + ' was replaced by ' + str(ale_interval[1]) + ' as the time step.\n'

    KratosMultiphysics.Logger.PrintInfo('MainStaticAngles', msg)

    # Adapt load restart settings
    if parameters_ale['solver_settings']['model_import_settings']['input_type'].GetString() == 'rest':
        folder_name = parameters_ale['solver_settings']['model_import_settings']['io_foldername'].GetString().split(path.sep)

        msg = 'Load restart settings adapted:'
        msg += '    - ' + folder_name[1] + ' was replaced by '
        msg += folder_ident + ' as the io_foldername\n'

        folder_name[1] = folder_ident
        parameters_ale['solver_settings']['model_import_settings']['io_foldername'].SetString(path.join(*folder_name))

        parameters_ale['solver_settings']['model_import_settings']['restart_load_file_label'].SetString(restart_time_label)

        msg += 'with restart from ' + restart_time_label + ' .'

        KratosMultiphysics.Logger.PrintInfo('MainStaticAngles', msg)

    # Adapting initial conditions
    for process in parameters_ale['processes']['initial_conditions_process_list']:
        # for a specific process update the z_rotation angle
        if process['python_module'].GetString() == "impose_mesh_motion_process":
            rot_z = process['Parameters']['z_rotation_angle'].GetDouble()

            process['Parameters']['z_rotation_angle'].SetDouble(float(static_case))

            msg = 'Rotation angle adapted ('+process['python_module'].GetString()+'):\n'
            msg += '    - ' + str(rot_z) + ' was replaced by ' + str(static_case) + ' as the z_rotation_angle.\n'
            KratosMultiphysics.Logger.PrintInfo('MainStaticAngles', msg)

    # Adapting auxiliary processes
    for process in parameters_ale['processes']['auxiliar_process_list']:
        if process['python_module'].GetString() == "single_mesh_temporal_output_process":

            file_name = process['Parameters']['file_settings']['file_name'].GetString().split(path.sep)

            msg = 'HDF5 file output adapted ('+process['python_module'].GetString()+'):\n'
            msg += '    - ' + file_name[1] + ' was replaced by ' + folder_ident
            msg += ' as the folder_name of files_settings.\n'

            file_name[1] = folder_ident
            process['Parameters']['file_settings']['file_name'].SetString(path.join(*file_name))

            # Updating time and step frequencies
            old_time_frequency = process['Parameters']['output_time_settings']['time_frequency'].GetDouble()
            old_step_frequency = process['Parameters']['output_time_settings']['step_frequency'].GetInt()

            new_step_frequency = 1
            new_time_frequency = 2*rampup_time

            process['Parameters']['output_time_settings']['time_frequency'].SetDouble(new_time_frequency)
            process['Parameters']['output_time_settings']['step_frequency'].SetInt(new_step_frequency)

            msg = 'HDF5 file output adapted ('+process['python_module'].GetString()+'):\n'
            msg += '    - ' + str(old_time_frequency) + ' was replaced by ' + str(new_time_frequency)
            msg += ' as the time_frequency.\n'
            msg += '    - ' + str(old_step_frequency) + ' was replaced by ' + str(new_step_frequency)
            msg += ' as the step_frequency.\n'

            KratosMultiphysics.Logger.PrintInfo('MainStaticAngles', msg)

        # Check that all processes have been adjusted
        else:
            msg = '[WARNING]: The process ' + process['python_module'].GetString() + ' does not match with any adjusted type.'
            KratosMultiphysics.Logger.PrintInfo('MainStaticAndles', msg)
            raise Exception('Unexpected auxiliary process. Make sure it is considered!')

    ###############################################################
    #
    # read and adapt CFD project parameters
    with open("ProjectParametersOpenMP_SA_CFD.json", 'r') as parameter_file:
        parameters_cfd = KratosMultiphysics.Parameters(parameter_file.read())

    print('\n---------- ADAPTING CFD PARAMETERS ----------\n')

    # Changing time parameters
    # Solver time added extra 10 time steps to solve for more than we ask output from

    old_start_time = parameters_cfd['problem_data']['start_time'].GetDouble()
    old_end_time = parameters_cfd['problem_data']['end_time'].GetDouble()
    old_delta_t = parameters_cfd['solver_settings']['time_stepping']['time_step'].GetDouble()
    new_start_time = cfd_interval[0]
    new_end_time = cfd_interval[1]

    parameters_cfd['problem_data']['start_time'].SetDouble(new_start_time)
    parameters_cfd['problem_data']['end_time'].SetDouble(new_end_time + 10*delta_t)
    parameters_cfd['solver_settings']['time_stepping']['time_step'].SetDouble(delta_t)

    msg = 'Time parameters adapted:\n'
    msg += '    - ' + str(old_start_time) + ' was replaced by ' + str(new_start_time) + ' as the start time.\n'
    msg += '    - ' + str(old_end_time) + ' was replaced by ' + str(new_end_time + 10*delta_t) + ' as the end time.\n'
    msg += '    - ' + str(old_delta_t) + ' was replaced by ' + str(delta_t) + ' as the time step.\n'
    KratosMultiphysics.Logger.PrintInfo('MainStaticAngles', msg)

    # Adapt save restart settings
    for process in parameters_cfd['output_processes']['my_processes']:
        if process['python_module'].GetString() == "save_restart_process":
            folder_name = process['Parameters']['io_foldername'].GetString().split(path.sep)

            msg = 'Save restart settings adapted ('+process['python_module'].GetString()+'):\n'
            msg += '    - ' + folder_name[1] + ' was replaced by '
            msg += folder_ident + 'as the io_foldername.\n'

            folder_name[1] = folder_ident
            process['Parameters']['io_foldername'].SetString(path.join(*folder_name))

            # Changing also the time frequency
            old_time_frequency = process['Parameters']['restart_save_frequency'].GetDouble()
            new_time_frequency = rampup_time

            process['Parameters']['restart_save_frequency'].SetDouble(new_time_frequency)

            msg += '    - ' + str(old_time_frequency) + ' was replaced by '
            msg += str(new_time_frequency) + ' as the time_frequency.\n'

            KratosMultiphysics.Logger.PrintInfo('MainStaticAngles', msg)

    # Adapting inlet (velocity parameters)
    for process in parameters_cfd['processes']['boundary_conditions_process_list']:
        if process['python_module'].GetString() == 'apply_inlet_process':

            old_interval = list(process['Parameters']['interval'].GetVector())
            old_modulus = process['Parameters']['modulus'].GetString()

            if 't' in process['Parameters']['modulus'].GetString():
                new_interval = rampup_interval
                new_modulus = str(U)+'*(t-'+str(rampup_interval[0])+')/'+str(rampup_time)
            else:
                # inlet boundary condition defined until 10 timesteps after simulation end
                new_interval = constant_condition_interval
                new_interval[1] += 20*delta_t
                new_modulus = str(U)

            process['Parameters']['interval'].SetVector(new_interval)
            process['Parameters']['modulus'].SetString(new_modulus)

            msg = 'Inlet process adapted ('+process['python_module'].GetString()+'):\n'
            msg += '    - ' + str(old_interval) + ' was replaced by ' + str(new_interval) + ' as the time interval.\n'
            msg += '    - ' + old_modulus + ' was relaced by ' + new_modulus + ' as the wind speed modulus.\n'

            KratosMultiphysics.Logger.PrintInfo('MainStaticAngles', msg)

    # Adapting various auxiliary processes
    for process in parameters_cfd['processes']['auxiliar_process_list']:

        # ASCII output
        if process['python_module'].GetString() in ['line_output_process',
                                                    'multiple_points_output_process',
                                                    'compute_body_fitted_drag_process',
                                                    'set_mesh_motion_and_get_forces_process',
                                                    'compute_global_force_process',
                                                    'cfl_output_process']:

            folder_name = process['Parameters']['output_file_settings']['folder_name'].GetString().split(path.sep)

            msg = 'Ascii file output adapted ('+process['python_module'].GetString()+'):\n'
            msg += '    - ' + folder_name[1] + ' was replaced by ' + folder_ident
            msg += ' as the folder_name of output_files_settings.\n'

            folder_name[1] = folder_ident
            process['Parameters']['output_file_settings']['folder_name'].SetString(path.join(*folder_name))

            KratosMultiphysics.Logger.PrintInfo('MainStaticAngles', msg)

            # rotate also output coordinates for certain processes
            # as point/line/multiples points output uses element seach by default
            # and element refers to the current coordinates (rotated)
            # and not the reference one (not rotated)
            if process['python_module'].GetString() == 'line_output_process':
                # only for certain types of line output
                if (('pressure' in process['Parameters']['output_file_settings']['folder_name'].GetString()
                    and 'contour' in process['Parameters']['output_file_settings']['folder_name'].GetString())
                    or ('velocity' in process['Parameters']['output_file_settings']['folder_name'].GetString()
                    and 'deck' in process['Parameters']['output_file_settings']['folder_name'].GetString())):

                    print(process['Parameters']['output_file_settings']['folder_name'].GetString())
                    print('Old')
                    print(process['Parameters']['start_point'].GetVector())
                    print(process['Parameters']['end_point'].GetVector())

                    for label in ["start_point", "end_point"]:

                        center = [0.0, 0.0]
                        cur_coord = process['Parameters'][label].GetVector()

                        new_xy = ccw_rotate_point_around_z(
                            center, cur_coord[:-1], math.radians(static_case))
                        new_coord = [new_xy[0], new_xy[1], cur_coord[2]]
                        process['Parameters'][label].SetVector(new_coord)

                    print('New')
                    print(process['Parameters']['start_point'].GetVector())
                    print(process['Parameters']['end_point'].GetVector())
                    #wait = input('check1...')

            elif process['python_module'].GetString() == 'multiple_points_output_process':
                # only for certain types of multiple points output
                if (('pressure' in process['Parameters']['output_file_settings']['folder_name'].GetString()
                    and 'coherence' in process['Parameters']['output_file_settings']['folder_name'].GetString())
                    or ('velocity' in process['Parameters']['output_file_settings']['folder_name'].GetString()
                    and 'selected' in process['Parameters']['output_file_settings']['folder_name'].GetString())):

                    positions = process['Parameters']['positions'].GetMatrix()

                    print(process['Parameters']['output_file_settings']['folder_name'].GetString())
                    print('Old')
                    print(process['Parameters']['positions'].GetMatrix())

                    num_points = positions.Size1()

                    center = [0.0, 0.0]
                    new_positions = KratosMultiphysics.Matrix(num_points, 3)
                    for k in range(num_points):

                        cur_coord = [0.0, 0.0, 0.0]
                        for j in range(3):
                            cur_coord[j] = positions[k,j]

                        new_xy = ccw_rotate_point_around_z(
                            center, cur_coord[:-1], math.radians(static_case))
                        new_coord = [new_xy[0], new_xy[1], cur_coord[2]]

                        new_positions[k,0] = new_coord[0]
                        new_positions[k,1] = new_coord[1]
                        new_positions[k,2] = new_coord[2]

                    process['Parameters']['positions'].SetMatrix(new_positions)

                    print('New')
                    print(process['Parameters']['positions'].GetMatrix())
                    #wait = input('check2...')

        # HDF5 output
        elif process['python_module'].GetString() in ['single_mesh_temporal_output_process']:

            file_name = process['Parameters']['file_settings']['file_name'].GetString().split(path.sep)

            msg = 'HDF5 file output adapted ('+process['python_module'].GetString()+'):\n'
            msg += '    - ' + file_name[1] + ' was replaced by ' + folder_ident
            msg += ' as the folder_name of files_settings.\n'

            file_name[1] = folder_ident
            process['Parameters']['file_settings']['file_name'].SetString(path.join(*file_name))

            # Updating time and step frequencies
            old_time_frequency = process['Parameters']['output_time_settings']['time_frequency'].GetDouble()
            old_step_frequency = process['Parameters']['output_time_settings']['step_frequency'].GetInt()

            # pressure output on deck
            if ('PRESSURE' == process['Parameters']['nodal_solution_step_data_settings']['list_of_variables'][0].GetString()
                    and process['Parameters']['nodal_solution_step_data_settings']['list_of_variables'].size() == 1 ):
                new_step_frequency = 1
            # all other output
            else:
                new_step_frequency = int(total_time/10/delta_t)

            new_time_frequency = delta_t*new_step_frequency

            msg += '    - ' + str(old_time_frequency) + ' was replaced by ' + str(new_time_frequency)
            msg += ' as the time_frequency.\n'
            msg += '    - ' + str(old_step_frequency) + ' was replaced by ' + str(new_step_frequency)
            msg += ' as the step_frequency.\n'

            process['Parameters']['output_time_settings']['time_frequency'].SetDouble(new_time_frequency)
            process['Parameters']['output_time_settings']['step_frequency'].SetInt(new_step_frequency)

            KratosMultiphysics.Logger.PrintInfo('MainStaticAngles', msg)

        # Statistics starting point
        elif process['python_module'].GetString() in ['temporal_statistics_process']:

            old_start_point = process['Parameters']['statistics_start_point_control_value'].GetDouble()
            new_start_time = converged_flow_interval[0]

            process['Parameters']['statistics_start_point_control_value'].SetDouble(new_start_time)

            msg = 'Adapting statistics process ('+process['python_module'].GetString()+'):\n'
            msg += '    - Replacing ' + str(old_start_point) + ' with ' + str(new_start_time) + ' as the statistics startig point.\n'
            KratosMultiphysics.Logger.PrintInfo('MainStaticAngles', msg)

        # Check that all processes have been adjusted
        else:
            msg = '[WARNING]: The process ' + process['python_module'].GetString() + ' does not match with any adjusted type.'
            KratosMultiphysics.Logger.PrintInfo('MainStaticAndles', msg)
            raise Exception('Unexpected auxiliary process. Make sure it is considered!')

    # needed as CFD used the model part from ALE and NOT the original MDPA
    process = parameters_cfd['solver_settings']['model_import_settings']
    process.RemoveValue('input_filename')
    process['input_type'].SetString("use_input_model_part")

    print('\n-------------- SETUP COMPLETED --------------\n')

    ###############################################################
    #
    # run first ALE, afterwards CFD
    model = KratosMultiphysics.Model()
    list_of_analyses = [
        MeshMovingAnalysis(model, parameters_ale),
        FluidDynamicsAnalysisWithFlush(model, parameters_cfd)
    ]

    for analysis in list_of_analyses:
        analysis.Run()
