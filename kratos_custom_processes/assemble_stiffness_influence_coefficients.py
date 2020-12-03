import os
import numpy as np

sub_dir = 'unsymmetric_y_1_floor.gid\\'
sub_dir_1 = 'unsymmetric_y\\'
src_dir = sub_dir + 'reaction_forces\\'

reaction_forces_files = os.listdir(src_dir)[9:] #if ground BC ground is in the folder
number_of_files = len(reaction_forces_files)
floors = int(number_of_files/9)
dofs = 3

#remove values close to zero 
def remove_small_values(matrix, tol = 1e-1):
    for i, row in enumerate(matrix):
        for j, val in enumerate(row):
            if abs(val) < tol:
                matrix[i][j] = 0.0

glob_K = np.zeros((floors * dofs,
                   floors * dofs))

begin_of_files = np.linspace(0,len(reaction_forces_files),dofs*dofs)
for j in range(floors):
    for i in range(floors):
        K_comb =[]
        for k in range(dofs):
            #print(reaction_forces_files[j*dofs*dofs + i*dofs + k])
            data = np.loadtxt(os.path.join(src_dir, reaction_forces_files[j*dofs*dofs + i*dofs + k]))
            K_comb.append([data[1], data[2], data[-1]])
        K_comb = np.transpose(np.asarray(K_comb))
        glob_K[i*dofs:i*dofs+dofs, dofs*j:dofs*j+dofs] += K_comb

remove_small_values(glob_K)

np.savetxt(sub_dir + "glob_K.csv", glob_K, delimiter= ' ')
print (glob_K)
