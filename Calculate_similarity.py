from quippy.potential import Potential
from quippy.descriptors import Descriptor
import numpy as np

# Import necessary libraries

# Define a function to calculate similarity between two sets of descriptors
def calculate_similarity(A, B):
    # Calculate the squared norms of descriptors in B
    B_norms = np.sum(B ** 2, axis=1)  
    # Compute the dot product between B and A
    AB = np.dot(B, A.T)  
    # Calculate the squared norms of descriptors in A
    A_norms = np.sum(A ** 2, axis=1)  
    # Compute the squared distances between descriptors
    distances_squared = B_norms[:, None] + A_norms - 2 * AB  
    # Set any negative distances to zero
    distances_squared[distances_squared < 0] = 0  
    # Calculate the distances by taking the square root
    distances = np.sqrt(distances_squared)  
    # Return the minimum distance for each descriptor in B
    return distances.min(axis=1)  

# Define a function to calculate the maximum and average similarities
def calculate_max_and_avg(similarities):
    # Calculate the maximum similarity value
    max_similarity = similarities.max()  
    # Calculate the average similarity value
    avg_similarity = similarities.mean()  
    # Return the maximum and average similarities
    return max_similarity, avg_similarity  
    
    
# Create a descriptor object using the "s_soap" descriptor
desc = Descriptor("s_soap")
# Compute the descriptor for the candidate structure
d1 = desc.calc(candidate_structure)
# Compute the descriptor for the training structure
d2 = desc.calc(training_structure)
# Calculate the similarities between the descriptors
similarities = calculate_similarity(d2["data"], d1["data"])
# Calculate the maximum and average similarities
Dmax, Davg = calculate_max_and_avg(similarities)