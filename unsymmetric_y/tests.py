# from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# import KratosMultiphysics
# from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
import os
import shutil
# with open("ProjectParametersCustom_JZ.json",'r') as parameter_file:
#     parameters = KratosMultiphysics.Parameters(parameter_file.read())
def order_files_in_folder(current_subdir_of_Kratos_stiffnessInfluence, dest_folder_name, file_suffix):
    subdir = current_subdir_of_Kratos_stiffnessInfluence
    dest_dir = os.path.join(subdir,dest_folder_name)
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
    for filename in os.listdir(subdir):
        if filename.endswith(file_suffix):
            shutil.move(os.path.join(subdir, filename), os.path.join(dest_dir, filename))

def move_and_rename_vtk_output (dest_folder_name, src_folder):
    subdir = 'unsymmetric_y'
    dest_folder_name = os.path.join(subdir, dest_folder_name)
    src_folder = os.path.join(subdir, src_folder)
    if not os.path.exists(dest_folder_name):
        os.makedirs(dest_folder_name)
    for filename in os.listdir(src_folder):
        #print(filename)
        os.rename(os.path.join(src_folder, filename), os.path.join(dest_folder_name, 'Structure_'+src_folder[-7:]+'.vtk'))
        os.rmdir(src_folder)

for dirs, subdirs, files in os.walk('unsymmetric_y'):
    for subdir in subdirs:
        if subdir.startswith('vtk_output_comb'):
            print (subdir)