# Importing the Kratos Library
import KratosMultiphysics

# other imports
from time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility


def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string")

    return ComputeIntegralValueProcess(model, settings["Parameters"])


class ComputeIntegralValueProcess(KratosMultiphysics.Process):
    """
    Auxiliary base class to output total flow forces
    over obstacles in fluid dynamics problems.
    A derived class needs to be implemented to be able to use
    this functionality, as calling the base class alone is not enough.
    """

    def __init__(self, model, params):
        """
        Auxiliary class to output total flow forces over obstacles
        in fluid dynamics problems for a body fitted model part.
        """
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"           : "",
                "interval"                  : [0.0, 1e30],
                "print_to_screen"      : false,
                "variable_name"        : "",
                "print_format"              : ".8f",
                "write_output_file"    : true,
                "output_file_settings": {}
            }
            """)

        self.kratos_vars = {}  # dict storing name-KratosVars,
        # hopefully faster than accessing KratosComponents all the time

        self.params = params
        # Detect "End" as a tag and replace it by a large number
        if(self.params.Has("interval")):
            if(self.params["interval"][1].IsString()):
                if(self.params["interval"][1].GetString() == "End"):
                    self.params["interval"][1].SetDouble(1e30)
                else:
                    raise Exception("The second value of interval can be \"End\" or a number, interval currently:" +
                                    self.params["interval"].PrettyPrintJsonString())

        self.params.ValidateAndAssignDefaults(default_settings)

        self.format = self.params["print_format"].GetString()

        # getting the ModelPart from the Model
        self.model_part_name = self.params["model_part_name"].GetString()
        if self.model_part_name == "":
            raise Exception('No "model_part_name" was specified!')
        else:
            self.model_part = model[self.model_part_name]

        self.var_name= self.params["variable_name"].GetString()
        if self.var_name == "":
            raise Exception('No "var_name" was specified!')
        else:
            self.var= self.__GetKratosVariable(self.var_name)

    def ExecuteInitialize(self):

        self.interval= KratosMultiphysics.Vector(2)
        self.interval[0]= self.params["interval"][0].GetDouble()
        self.interval[1]= self.params["interval"][1].GetDouble()
        self.print_to_screen= self.params["print_to_screen"].GetBool()
        self.write_output_file= self.params["write_output_file"].GetBool()

        if (self.model_part.GetCommunicator().MyPID() == 0):
            if (self.write_output_file):
                file_handler_params= KratosMultiphysics.Parameters(
                    self.params["output_file_settings"])

                file_header=self._GetFileHeader()
                self.output_file=TimeBasedAsciiFileWriterUtility(self.model_part,
                                                                   file_handler_params, file_header).file

    def ExecuteFinalizeSolutionStep(self):

        current_time=self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if((current_time >= self.interval[0]) and (current_time < self.interval[1])):
            # Compute the force
            integral_value=self._GetCorrespondingIntegralValue()

            # Write the value components
            if (self.model_part.GetCommunicator().MyPID() == 0):
                if (self.print_to_screen):
                    result_msg='Integral value for model part ' + \
                        self.model_part_name + ' and variable ' + self.var_name + '\n'
                    result_msg += str(current_time) + " x-comp: " + format(integral_value[0], self.format) + " y-comp: " + format(
                        integral_value[1], self.format) + " z-comp: " + format(integral_value[2], self.format)
                    self._PrintToScreen(result_msg)

                # not formatting time in order to not lead to problems with time recognition
                # in the file writer when restarting
                if (self.write_output_file):
                    self.output_file.write(str(current_time)+" "+format(integral_value[0], self.format)+" "+format(
                        integral_value[1], self.format)+" "+format(integral_value[2], self.format)+"\n")

    def ExecuteFinalize(self):
        if (self.model_part.GetCommunicator().MyPID() == 0):
            self.output_file.close()

    def _GetFileHeader(self):
        header = '# Integral value for model part ' + self.model_part_name + ' and variable ' + self.var_name + '\n'
        header += '# Time VectorVal[0] VectorVal[1] VectorVal[2]\n'
        return header

    def _PrintToScreen(self, result_msg):
        KratosMultiphysics.Logger.PrintInfo(
            "ComputeIntegralValueProcess", "INTEGRAL VALUE RESULTS:")
        KratosMultiphysics.Logger.PrintInfo(
            "ComputeIntegralValueProcess", "Current time: " + result_msg)

    def _GetCorrespondingIntegralValue(self):
        # TODO: for now no type-check, should be checked if vector or not...
        '''
        SOMETHING LIKE THIS:


            if type(kratos_var) == KratosMultiphysics.DoubleVariable or type(kratos_var) == KratosMultiphysics.Array1DComponentVariable:
                SetData(model_part, kratos_var, data_array)
            elif type(kratos_var) == KratosMultiphysics.Array1DVariable3:
                domain_size = model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
                if not domain_size in [1,2,3]:
                    raise Exception("DOMAIN_SIZE has to be 1, 2 or 3!")
                num_nodes = NumberOfNodes(model_part)
                if data_array.size != num_nodes*domain_size:
                    raise Exception("Size of data does not match number of nodes x domain size!")
                ext = ["_X", "_Y", "_Z"]
                for i in range(domain_size):
                    component_var = self.__GetKratosVariable(kratos_var.Name()+ext[i])
                    range_begin = i*num_nodes
                    range_end = (i+1)*num_nodes
                    SetData(model_part, component_var, data_array[range_begin:range_end])
            else:
                err_msg  = 'Type of variable "' + kratos_var.Name() + '" is not valid\n'
                err_msg += 'It can only be double, component or array3d!'
                raise Exception(err_msg)

        '''
        vector_val = [0.0, 0.0, 0.0]

        for node in self.model_part.Nodes:
            nodal_result = node.GetSolutionStepValue(self.var, 0)

            vector_val[0] += nodal_result[0]
            vector_val[1] += nodal_result[1]
            vector_val[2] += nodal_result[2]

        return vector_val

    def __GetKratosVariable(self, var_name):
        if not var_name in self.kratos_vars:
            self.kratos_vars[var_name] = KratosMultiphysics.KratosGlobals.GetVariable(
                var_name)

        return self.kratos_vars[var_name]
