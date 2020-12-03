# Importing the Kratos Library
import KratosMultiphysics

# other imports
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility


def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return ComputeCustomBaseReactionProcess(Model, settings["Parameters"])


class ComputeCustomBaseReactionProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"           : "please_specify_model_part_name",
                "interval"                  : [0.0, 1e30],
                "reference_point"           : [0.0, 0.0, 0.0],
                "print_drag_to_screen"      : false,
                "print_format"              : ".8f",
                "write_drag_output_file"    : true,
                "output_file_settings": {}
            }
            """)

        self.settings = settings

        # Detect "End" as a tag and replace it by a large number
        if(self.settings.Has("interval")):
            if(self.settings["interval"][1].IsString()):
                if(self.settings["interval"][1].GetString() == "End"):
                    self.settings["interval"][1].SetDouble(1e30)
                else:
                    raise Exception("The second value of interval can be \"End\" or a number, interval currently:"+self.settings["interval"].PrettyPrintJsonString())

        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model[self.settings["model_part_name"].GetString()]
        self.interval = KratosMultiphysics.Vector(2)
        self.interval[0] = self.settings["interval"][0].GetDouble()
        self.interval[1] = self.settings["interval"][1].GetDouble()
        self.print_drag_to_screen = self.settings["print_drag_to_screen"].GetBool()
        self.write_drag_output_file = self.settings["write_drag_output_file"].GetBool()

        self.format = self.settings["print_format"].GetString()

        # PMT: added reference point for moment calculation
        self.reference_x = self.settings["reference_point"][0].GetDouble()
        self.reference_y = self.settings["reference_point"][1].GetDouble()
        self.reference_z = self.settings["reference_point"][2].GetDouble()

        if (self.model_part.GetCommunicator().MyPID() == 0):
            if (self.write_drag_output_file):

                file_handler_params= KratosMultiphysics.Parameters(
                    settings["output_file_settings"])

                file_header=self._GetFileHeader()
                self.output_file=TimeBasedAsciiFileWriterUtility(self.model_part,
                                                                   file_handler_params, file_header).file


    def ExecuteFinalizeSolutionStep(self):

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if((current_time >= self.interval[0]) and  (current_time < self.interval[1])):

            # Note that MPI communication is done within VariableUtils().SumHistoricalNodeVectorVariable()
            #reaction_vector = KratosMultiphysics.VariableUtils().SumHistoricalNodeVectorVariable(KratosMultiphysics.REACTION, self.model_part, 0)

            # PMT: TODO: only checked for OpenMP, update for MPI
            fx = 0.0
            fy = 0.0
            fz = 0.0

            mx = 0.0
            my = 0.0
            mz = 0.0       
            
            for node in self.model_part.Nodes:
                reaction =        node.GetSolutionStepValue(KratosMultiphysics.REACTION, 0)
                moment_reaction = node.GetSolutionStepValue(KratosMultiphysics.REACTION_MOMENT, 0)

                # PMT: NOTE: sign is flipped to go from reaction to action
                # NOTE flipped it back
                fx += (1) * reaction[0]
                fy += (1) * reaction[1]
                fz += (1) * reaction[2]

                #print (node.id, reaction[0])

                x = node.X - self.reference_x
                y = node.Y - self.reference_y
                z = node.Z - self.reference_z
                mx += y * (1) * reaction[2] - z * (1) * reaction[1] + (1) * moment_reaction[0]
                my += z * (1) * reaction[0] - x * (1) * reaction[2] + (1) * moment_reaction[1]
                mz += x * (1) * reaction[1] - y * (1) * reaction[0] + (1) * moment_reaction[2]
                
            if (self.model_part.GetCommunicator().MyPID() == 0):

                if (self.print_drag_to_screen):
                    print("CUSTOM BASE REACTION RESULTS:")
                    print("Current time: " + str(current_time))
                    print("Forces:" + " Fx: " + str(fx) + " Fy: " + str(fy) + " Fz: " + str(fz))
                    print("Moments:" + " Mx: " + str(mx) + " My: " + str(my) + " Mz: " + str(mz))

                if (self.write_drag_output_file):
                    #     with open(self.drag_filename, 'a') as file:
                    #         output_str = str(current_time)
                    #         output_str += "   " + str(fx) + "   " + str(fy) + "   " + str(fz)
                    #         output_str += "   " + str(mx) + "   " + str(my) + "   " + str(mz) + "\n"
                    #         file.write(output_str)
                    #         file.close()


                    # if (self.write_output_file):
                    # output_str = str(current_time)
                    # output_str += "   " + format(fx, self.format) + "   " + format(fy, self.format) + "   " + format(fz, self.format)
                    # output_str += "   " + format(mx, self.format) + "   " + format(my, self.format) + "   " + format(mz, self.format) + "\n"
                    
                    output_str = str(current_time)
                    output_str += "   " + str(fx) + "   " + str(fy) + "   " + str(fz)
                    output_str += "   " + str(mx) + "   " + str(my) + "   " + str(mz) + "\n"

                    self.output_file.write(output_str)

                    # self.output_file.write(str(current_time)+" "+format(integral_value[0], self.format)+" "+format(
                    #     integral_value[1], self.format)+" "+format(integral_value[2], self.format)+"\n")

    def _GetFileHeader(self):
        header = '# Integral value for model part ' + self.settings["model_part_name"].GetString() + ' ' + 'at reference point X = ' + str(self.reference_x) + ', Y = ' + str(self.reference_y)+ ', Z = ' + str(self.reference_z) + '\n'
        header += '# Time Fx: Fy: Fz: Mx: My: Mz:\n'
        return header

    def ExecuteFinalize(self):
        if (self.model_part.GetCommunicator().MyPID() == 0):
            self.output_file.close()