from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7


import sys
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
import math
import shutil
import os
"""
For user-scripting it is intended that a new class is derived
from StructuralMechanicsAnalysis to do modifications
"""
def continue_run():
    stop = input('Want to continue? (y/n) ')
    if stop == 'n':
        sys.exit()

def order_files_in_folder(dest_folder_name, file_suffix):
    dest_dir = dest_folder_name
    if not os.path.exists(dest_folder_name):
        os.makedirs(dest_dir)
    for filename in os.listdir():
        if filename.endswith(file_suffix):
            shutil.move(filename, os.path.join(dest_dir, filename))

def move_and_rename_vtk_output (dest_folder_name):
    if not os.path.exists(dest_folder_name):
        os.makedirs(dest_folder_name)
    for dirs, subdirs, files in os.walk(os.getcwd()):
        for subdir in subdirs:
            if subdir.startswith('vtk_output_comb'):
                for filename in os.listdir(subdir):
                    src = os.path.join(subdir, filename)
                    dest = os.path.join(dest_folder_name, 'Structure_' + subdir[-7:] + '.vtk')
                    os.rename(src, dest)
                    os.rmdir(subdir)

if __name__ == "__main__":

    '''
    the base combinations are applied one after the other
    x displacement 1 m
    y displacement 1 m 
    z rotation 1 rad
    '''
    base_combinations = [[1.0,0.0],[0.0,1.0],1.0,[0.0,0.0]] # x = 1, y= 1, theta = 1 °
    fixed = [[1,2],[-1,1],[-1,2]]
    combinations = {}

    for floor in range(3):
        for dof in range(3):
            if floor == 0:
                fixed_1, fixed_2 = 1,2
            elif floor == 1:
                fixed_1, fixed_2 = -1,1
            elif floor == 2:
                fixed_1, fixed_2 = -1,-2
            if dof == 2:
                combinations['comb_'+str(floor)+str(dof)] = {
                    "Structure.DISPLACEMENT_plate_"+str(floor):{"dxdy": [0.0, 0.0], 'dtheta':base_combinations[dof]},
                    "Structure.DISPLACEMENT_plate_"+str(floor+fixed_1):{"dxdy": [0.0, 0.0], 'dtheta':0.0},
                    "Structure.DISPLACEMENT_plate_"+str(floor+fixed_2):{"dxdy": [0.0, 0.0], 'dtheta':0.0}
                    }        
            else:
                combinations['comb_'+str(floor)+str(dof)] = {
                    "Structure.DISPLACEMENT_plate_"+str(floor):{'dxdy':base_combinations[dof], 'dtheta':0.0},
                    "Structure.DISPLACEMENT_plate_"+str(floor+fixed_1):{'dxdy':base_combinations[-1], 'dtheta':0.0},
                    "Structure.DISPLACEMENT_plate_"+str(floor+fixed_2):{'dxdy':base_combinations[-1], 'dtheta':0.0}
                    }
                
    groups_list = ["Structure.DISPLACEMENT_plate_"+str(floor) for floor in range(3)]

    #print('The combinations are: \n' , combinations)
    
    # run a simulation for each combination 
    for my_comb in combinations:
        
        with open("ProjectParametersCustom_1.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        # only mesh moving process is used 
        for process in parameters['processes']['constraints_process_list']:
            
            # displacements x and y: dof 0 and 1
            if process['python_module'].GetString() == 'impose_mesh_motion_process':
                if my_comb[-1] != '2': # relative displacement x and y 
                    for my_group in combinations[my_comb]:
                        if my_group == process['Parameters']['model_part_name'].GetString():
                            process['Parameters']['dxdy'].SetVector(combinations[my_comb][my_group]['dxdy'])
                # dof 2 = rotation around z            
                elif my_comb[-1] == '2': 
                    for my_group in combinations[my_comb]:
                        if my_group == process['Parameters']['model_part_name'].GetString():
                            process['Parameters']['dtheta'].SetDouble(combinations[my_comb][my_group]['dtheta'])

        # for gid output adjust file name
        for process in parameters['output_processes']['gid_output']:
            if process['python_module'].GetString() == 'gid_output_process':
            
                process['Parameters']['output_name'].SetString(process['Parameters']['output_name'].GetString() + "_" + my_comb)

                #wait = input('check gid output naming...')
            else:
                Exception("Not a valid gid output process")

        # for vtk output adjust folder name
        for process in parameters['output_processes']['vtk_output']:
            if process['python_module'].GetString() == 'vtk_output_process':
                src_folder = process['Parameters']['folder_name'].GetString() + "_" + my_comb
                process['Parameters']['folder_name'].SetString(src_folder)
                #process['Parameters']['model_part_name'].SetString(process['Parameters']['model_part_name'].GetString() + "_" + my_comb)

                #wait = input('check vtk output naming...')
            else:
                Exception("Not a valid vtk output process")
        
        # # Add custom base reaction processes for each group and maybe for each node one up
        
        base_reaction_process_list = []
        for group in groups_list:
            base_reaction_process_list.append(
                {
                "python_module" : "compute_custom_base_reaction_process",
                "process_name"  : "ComputeCustomBaseReactionProcess",
                "Parameters"    : {
                    "model_part_name"        : group,
                    "write_drag_output_file" : True,
                    "print_drag_to_screen"   : False,
                    "reference_point"        : [0.0, 0.0, 0.0],
                    "interval"               : [0.0,"End"],
                    "output_file_settings": {
                        "file_name"   : "base_reaction" + "_" + group[-7:]}
                }
            }
            )
        # parameters['processes']['list_other_processes']:base_reaction_process_list
        
        # for the base reaction output adjust file name
        for process in parameters['processes']['list_other_processes']:
            if process['python_module'].GetString() == 'compute_custom_base_reaction_process':
                # group = process['Parameters']['model_part_name'][-7:]
                #print("base reaction output process - current name: ", process['Parameters']['output_file_settings']['file_name'].GetString())
                process['Parameters']['output_file_settings']['file_name'].SetString(
                    process['Parameters']['output_file_settings']['file_name'].GetString() + "_" + my_comb)
                
                #print("base reaction output process - new name: ", process['Parameters']['output_file_settings']['file_name'].GetString())

                #wait = input('check base reacion output naming...')

        # for each combination run Kratos
        model = KratosMultiphysics.Model()
        simulation = StructuralMechanicsAnalysis(model,parameters)
        simulation.Run()

    order_files_in_folder('gid_post_bins', '.bin')
    order_files_in_folder('reaction_forces', '.dat')
    move_and_rename_vtk_output('vtk_output_all')