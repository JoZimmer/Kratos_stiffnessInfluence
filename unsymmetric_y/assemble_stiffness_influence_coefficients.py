import os
import numpy as np

print(os.getcwd())
# src_dir = 'unsymmetric_y\\reaction_forces\\'
# reaction_forces_files = os.listdir(src_dir)[9:] #if ground BC is in the folder

# #remove values close to zero 
# def remove_small_values(matrix, tol = 1e-1):
#     for i, row in enumerate(matrix):
#         for j, val in enumerate(row):
#             if abs(val) < tol:
#                 matrix[i][j] = 0.0

# K = []
# for k in range(0,27,9):
#     K_plate = []
#     for j in range(3):
#         K_comb = []
#         for i in range(3):
#             data = np.loadtxt(os.path.join(src_dir, reaction_forces_files[k + 3*j + i]))
#             K_comb.append([data[1], data[2], data[-1]])
#         K_comb = np.transpose(np.asarray(K_comb))
#         remove_small_values(K_comb)    
#         K_plate.append(K_comb)
#     K_plate = np.concatenate((K_plate[0],K_plate[1],K_plate[2]), axis = 0)
#     K.append(K_plate)
# K = np.concatenate((K[0],K[1],K[2]), axis = 1)

# np.savetxt("unsymmetric_y\\K.csv", K, delimiter= ' ')
# print (K)
