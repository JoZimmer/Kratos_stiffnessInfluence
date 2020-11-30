import KratosMultiphysics
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility
import sys
import math


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


def ccw_rotate_vect_comp_around_z(vect, angle):
    '''
    Counter-clockwise roration around z-axis in a right-handed coordinate system
    '''
    rot_matr_z = [[math.cos(angle), math.sin(angle), 0],
                  [-math.sin(angle), math.cos(angle), 0],
                  [0, 0, 1]]

    vect_rot = [0.0, 0.0, 0.0]

    for i in range(3):
        vect_rot[i] = sum([a*b for a, b in zip(rot_matr_z[i], vect)])

    return vect_rot


def Factory(params, Model):
    if(type(params) != KratosMultiphysics.Parameters):
        raise Exception(
            'expected input shall be a Parameters object, encapsulating a json string')
    return SetMeshMotionAndGetForcesProcess(Model, params["Parameters"])


class SetMeshMotionAndGetForcesProcess(KratosMultiphysics.Process):
    '''
    Computes the flow- and body-attached forces
    for a model part with body-fitted mesh
    split into a number of intervals
    thus the naming LevelForces

    Takes as input a CCW positive rotation (in degrees) around
    axis z for the body-attached forces
    '''

    def __init__(self, Model, params):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"       : "",
                "interval"              : [0.0, 1e30],
                "rampup_time"           : 0.0,
                "reference_point"       : [0.0,0.0,0.0],
                "z_rotation_angle"      : 0.0,
                "imposed_motion":{
                    "pitch": {"amplitude": [0.1], "frequency" : [1.2]},
                    "heave": {"amplitude": [0.02, 0.03], "frequency" : [0.9, 1.1]}
                },
                "print_to_screen"       : false,
                "print_format"          : ".8f",
                "write_output_file"     : true,
                "output_file_settings"  : {}
            }
            """)

        # Detect 'End' as a tag and replace it by a large number
        if(params.Has('interval')):
            if(params['interval'][1].IsString()):
                if(params['interval'][1].GetString() == 'End'):
                    params['interval'][1].SetDouble(1e30)
                else:
                    raise Exception('The second value of interval can be \'End\' or a number, interval currently:' +
                                    params['interval'].PrettyPrintJsonString())

        params.ValidateAndAssignDefaults(default_settings)

        # modal part params
        self.model_part_name = params['model_part_name'].GetString()
        self.model_part = Model[self.model_part_name]
        self.interval = params["interval"].GetVector()
        self.print_to_screen = params['print_to_screen'].GetBool()
        self.write_output_file = params['write_output_file'].GetBool()
        self.format = params["print_format"].GetString()
        self.rampup_time = params['rampup_time'].GetDouble()

        # added reference point for moment calculation
        reference = params['reference_point'].GetVector()
        if reference.Size() != 3:
            raise Exception(
                'The reference point position has to be provided with 3 coordinates!')
        self.reference_center = reference
        self.updated_center = reference

        # user inpput expected in degrees, here changing to radians
        self.reference_theta = math.radians(
            params['z_rotation_angle'].GetDouble())
        self.updated_theta = self.reference_theta

        self.motion_increment = {"pitch": 0.0, "heave": 0.0}
        self.prescribed_motion = {
            "pitch": {"amplitude": params['imposed_motion']['pitch']['amplitude'].GetVector(),
                      "frequency": params['imposed_motion']['pitch']['frequency'].GetVector()},
            "heave": {"amplitude": params['imposed_motion']['heave']['amplitude'].GetVector(),
                      "frequency": params['imposed_motion']['heave']['frequency'].GetVector()}}

        self.prescribed_motion['pitch']['amplitude'] = [math.radians(a) for a in self.prescribed_motion['pitch']['amplitude']]

        if (self.model_part.GetCommunicator().MyPID() == 0):
            if (self.write_output_file):

                # default name as fallback
                output_file_name = params["model_part_name"].GetString()

                file_handler_params = KratosMultiphysics.Parameters(
                    params["output_file_settings"])

                self.output_file = {}
                for case in ['motion', 'force']:
                    case_file_handler_params = file_handler_params.Clone()

                    if file_handler_params.Has("file_name"):
                        output_file_name = file_handler_params["file_name"].GetString(
                        )
                        warn_msg = 'Unexpected user-specified entry found in "output_file_settings": {"file_name": '
                        warn_msg += '"' + \
                            file_handler_params["file_name"].GetString(
                            ) + '"}\n'
                        warn_msg += 'Using this specififed file name instead of the default "' + \
                            output_file_name + '"'
                        KratosMultiphysics.Logger.PrintWarning(
                            "SetMeshMotionAndGetForcesProcess", warn_msg)
                    else:
                        case_file_handler_params.AddEmptyValue("file_name")

                    case_file_handler_params["file_name"].SetString(
                        output_file_name + '_' + case + '.dat')

                    file_header = self._GetFileHeader(case)
                    self.output_file[case] = TimeBasedAsciiFileWriterUtility(self.model_part,
                                                                             case_file_handler_params, file_header).file

    def ExecuteInitializeSolutionStep(self):
        self._UpdateSolutionStepIncrement()

        for node in self.model_part.Nodes:
            # using X0 and Y0 - initial coordinates - as deformation prescribed
            # with respect to this undeformed state
            cur_x0y0 = [node.X0, node.Y0]
            new_x0y0 = ccw_rotate_point_around_z(
                self.reference_center[:-1], cur_x0y0, self.updated_theta)

            # calculate and set increment
            # due to pitch
            [dx, dy] = [a-b for a, b in zip(new_x0y0, cur_x0y0)]
            # due to heave
            dy += self.motion_increment['heave']

            node.SetSolutionStepValue(
                KratosMultiphysics.MESH_DISPLACEMENT_X, dx)
            node.SetSolutionStepValue(
                KratosMultiphysics.MESH_DISPLACEMENT_Y, dy)

    def _UpdateSolutionStepIncrement(self):
        if self.model_part.ProcessInfo[KratosMultiphysics.TIME] >= self.rampup_time:

            for key in self.motion_increment.keys():
                # using list syntax
                # such that iteratively superposition
                # of prescribed amplitudes and frequencies is possible
                # motion increment computed with respect to reference state
                # val = sum_of_contributions(a_i*sin(2*pi*f*(t-t_ramp)))
                self.motion_increment[key] = sum([a*math.sin(2*math.pi*f*(self.model_part.ProcessInfo[KratosMultiphysics.TIME]-self.rampup_time))
                                                  for a, f in zip(self.prescribed_motion[key]['amplitude'],
                                                                  self.prescribed_motion[key]['frequency'])])

        # heave -> needs updating for force calculation
        self.updated_center[1] = self.reference_center[1] + \
            self.motion_increment['heave']

        # pitch -> does not need updating for force calculation as moment does not depent on it
        # it is only important for motion
        self.updated_theta = self.reference_theta + \
            self.motion_increment['pitch']

    def ExecuteFinalizeSolutionStep(self):

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if((current_time >= self.interval[0]) and (current_time < self.interval[1])):
            ff, mf, fb, mb = self._EvaluateGlobalForces()

            if (self.model_part.GetCommunicator().MyPID() == 0):
                output_f = []
                output_f.extend(ff)
                output_f.extend(mf)
                output_f.extend(fb)
                output_f.extend(mb)

                output_f_vals = [format(val, self.format) for val in output_f]
                output_m_vals = [format(val, self.format) for val in [
                    0.0, self.motion_increment['heave'], self.motion_increment['pitch']]]
                # not formatting time in order to not lead to problems with time recognition
                # in the file writer when restarting
                output_f_vals.insert(0, str(current_time))
                output_m_vals.insert(0, str(current_time))

                res_f_labels = ['time: ',
                                'fx: ', 'fy: ', 'fz: ', 'mx: ', 'my: ', 'mz: ',
                                'fx\': ', 'fy\': ', 'fz\': ', 'mx\': ', 'my\': ', 'mz\': ']

                res_m_labels = ['time: ',
                                'dx: ', 'dy: ', 'dtheta: ']

                if (self.print_to_screen):
                    result_msg = []

                    result_msg.append('Force evaluation for model part ' +
                                      self.model_part_name + '\n')
                    result_msg[0] += ', '.join([a+b for a,
                                                b in zip(res_f_labels, output_f_vals)])

                    result_msg.append('Prescribed motion for model part ' +
                                      self.model_part_name + '\n')
                    result_msg[1] += ', '.join([a+b for a,
                                                b in zip(res_m_labels, output_m_vals)])

                    self._PrintToScreen(result_msg)
                    sys.stdout.flush()

                if (self.write_output_file):
                    self.output_file['force'].write(
                        ' '.join(output_f_vals) + '\n')
                    self.output_file['motion'].write(
                        ' '.join(output_m_vals) + '\n')

    def _EvaluateGlobalForces(self):
        # flow-attached forces: in x-y-z coordinate system
        ff = [0.0, 0.0, 0.0]
        mf = [0.0, 0.0, 0.0]

        for node in self.model_part.GetCommunicator().LocalMesh().Nodes:
            # sign is flipped to go from reaction to action -> force
            nodal_force = (-1) * \
                node.GetSolutionStepValue(KratosMultiphysics.REACTION, 0)

            # summing up nodal contributions to get resultant for model_part
            ff[0] += nodal_force[0]
            ff[1] += nodal_force[1]
            ff[2] += nodal_force[2]

            # using X, Y, Z as forceas are calculated
            # with respect to current (deformed) configuration
            x = node.X - self.updated_center[0]
            y = node.Y - self.updated_center[1]
            z = node.Z - self.updated_center[2]
            mf[0] += y * nodal_force[2] - z * nodal_force[1]
            mf[1] += z * nodal_force[0] - x * nodal_force[2]
            mf[2] += x * nodal_force[1] - y * nodal_force[0]

        # body-attached forces -> here only a rotation around z-axis
        # in x'-y'-z' coordinate system
        # of the summed-up forces and moment

        fb = ccw_rotate_vect_comp_around_z(ff, self.updated_theta)
        mb = ccw_rotate_vect_comp_around_z(mf, self.updated_theta)

        ff = self.model_part.GetCommunicator().GetDataCommunicator().SumDoubles(ff, 0)
        mf = self.model_part.GetCommunicator().GetDataCommunicator().SumDoubles(mf, 0)

        fb = self.model_part.GetCommunicator().GetDataCommunicator().SumDoubles(fb, 0)
        mb = self.model_part.GetCommunicator().GetDataCommunicator().SumDoubles(mb, 0)

        return ff, mf, fb, mb

    def _GetFileHeader(self, case):
        header = {}

        header['force'] = '# Forces for model part ' + \
            self.model_part_name + '\n'
        header['force'] += '# Time Fx Fy Fz Mx My Mz Fx\' Fy\' Fz\' Mx\' My\' Mz\'\n'

        header['motion'] = '# Motion for model part ' + \
            self.model_part_name + '\n'
        header['motion'] += '# Time dx dy dtheta\n'

        return header[case]

    def _PrintToScreen(self, result_msg):
        KratosMultiphysics.Logger.PrintInfo(
            'SetMeshMotionAndGetForcesProcess', 'Force - flow- and body-attached:')
        KratosMultiphysics.Logger.PrintInfo(
            'SetMeshMotionAndGetForcesProcess', 'Current time: ' + result_msg[0])
        KratosMultiphysics.Logger.PrintInfo(
            'SetMeshMotionAndGetForcesProcess', 'Motion:')
        KratosMultiphysics.Logger.PrintInfo(
            'SetMeshMotionAndGetForcesProcess', 'Current time: ' + result_msg[1])
