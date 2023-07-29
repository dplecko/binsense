
import numpy as np
import cpp_scaling_2d

# Generate a random matrix Sigma of size 4x4
Sigma = np.random.rand(4, 4)

# Apply the C++ function to the matrix
result = cpp_scaling_2d.cpp_scaling_2d(Sigma)
print(result)
