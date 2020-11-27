from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

"""
For user-scripting it is intended that a new class is derived
from StructuralMechanicsAnalysis to do modifications
"""

if __name__ == "__main__":

    combinations = {
        "comb1": {
            "Structure.DISPLACEMENT_Group0" : {
                "dxdy": [0.0, 0.0],
                "dtheta": 0.0,
                "x0y0" : [0.0,0.0]
            },
            "Structure.DISPLACEMENT_Group2" : {
                "dxdy": [1.2, 2.1],
                "dtheta": 10,
                "x0y0" : [0.0, 0.0]
            }
        },
        "comb2": {
            "Structure.DISPLACEMENT_Group0" : {
                "dxdy": [1.2, 2.1],
                "dtheta": 10,
                "x0y0" : [0.0, 0.0]
            },
            "Structure.DISPLACEMENT_Group2" : {
                "dxdy": [0.0, 0.0],
                "dtheta": 0.0,
                "x0y0" : [0.0,0.0]}
        }
    }

    for my_comb in combinations:

        with open("ProjectParametersCustom.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        # adapting the (mesh) motion process
        for process in parameters['processes']['constraints_process_list']:
            if process['python_module'].GetString() == 'impose_mesh_motion_process':
                for my_group in combinations[my_comb]:
                    if my_group == process['Parameters']['model_part_name'].GetString():
                        #print(combinations[my_comb][my_group]['dxdy'])
                        process['Parameters']['dxdy'].SetVector(combinations[my_comb][my_group]['dxdy'])
                        #print(combinations[my_comb][my_group]['dtheta'])
                        process['Parameters']['dtheta'].SetDouble(combinations[my_comb][my_group]['dtheta'])
                        #print(combinations[my_comb][my_group]['x0y0'])
                        process['Parameters']['x0y0'].SetVector(combinations[my_comb][my_group]['x0y0'])

                        #wait = input('check impose mesh motion process input...')

        # for gid output adjust file name
        for process in parameters['output_processes']['gid_output']:
            if process['python_module'].GetString() == 'gid_output_process':
                
                #print("gid output process - current name: ", process['Parameters']['output_name'].GetString())
                process['Parameters']['output_name'].SetString(process['Parameters']['output_name'].GetString() + "_" + my_comb)
                #print("gid output process - new name: ", process['Parameters']['output_name'].GetString())

                #wait = input('check gid output naming...')
            else:
                Exception("Not a valid gid output process")

        # for vtk output adjust folder name
        for process in parameters['output_processes']['vtk_output']:
            if process['python_module'].GetString() == 'vtk_output_process':
                
                #print("vtk output process - current name: ", process['Parameters']['folder_name'].GetString())
                process['Parameters']['folder_name'].SetString(process['Parameters']['folder_name'].GetString() + "_" + my_comb)
                #print("vtk output process - new name: ", process['Parameters']['folder_name'].GetString())

                #wait = input('check vtk output naming...')
            else:
                Exception("Not a valid vtk output process")

        # for the base reaction output adjust file name
        for process in parameters['processes']['list_other_processes']:
            if process['python_module'].GetString() == 'compute_custom_base_reaction_process':
                
                #print("base reaction output process - current name: ", process['Parameters']['output_file_settings']['file_name'].GetString())
                process['Parameters']['output_file_settings']['file_name'].SetString(process['Parameters']['output_file_settings']['file_name'].GetString() + "_" + my_comb)
                #print("base reaction output process - new name: ", process['Parameters']['output_file_settings']['file_name'].GetString())

                #wait = input('check base reacion output naming...')
    
        model = KratosMultiphysics.Model()
        simulation = StructuralMechanicsAnalysis(model,parameters)
        simulation.Run()
