import KratosMultiphysics
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


def Factory(params, Model):
    if(type(params) != KratosMultiphysics.Parameters):
        raise Exception(
            'expected input shall be a Parameters object, encapsulating a json string')
    return ImposeMeshMotionProcess(Model, params["Parameters"])


class ImposeMeshMotionProcess(KratosMultiphysics.Process):
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
                "z_rotation_angle"      : 0.0,
                "reference_point"       : [0.0,0.0,0.0],
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

        # added reference point for moment calculation
        reference = params['reference_point'].GetVector()
        if reference.Size() != 3:
            raise Exception(
                'The reference point position has to be provided with 3 coordinates!')
        self.reference_center = reference

        # user inpput expected in degrees, here changing to radians
        self.theta = math.radians(params['z_rotation_angle'].GetDouble())

        # alternative (if you used my suggestion):
        # self.theta = math.radians(self.model_part.GetRootModelPart[KratosMultiphysics.ROTATION])

    def ExecuteInitializeSolutionStep(self):
        for node in self.model_part.Nodes:
            # using X0 and Y0 - initial coordinates - as deformation prescribed
            # with respect to this undeformed state
            cur_x0y0 = [node.X0, node.Y0]
            new_x0y0 = ccw_rotate_point_around_z(
                self.reference_center[:-1], cur_x0y0, self.theta)

            # calculate and set increment
            # due to pitch
            [dx, dy] = [a-b for a, b in zip(new_x0y0, cur_x0y0)]

            node.SetSolutionStepValue(
                KratosMultiphysics.MESH_DISPLACEMENT_X, dx)
            node.SetSolutionStepValue(
                KratosMultiphysics.MESH_DISPLACEMENT_Y, dy)
