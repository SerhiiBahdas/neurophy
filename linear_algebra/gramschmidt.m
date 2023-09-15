function orthonormal_vectors = gramschmidt(vectors)

%GRAMSCHMIDT Orthogonalize a set of vectors using the Gram-Schmidt process.
%
%   orthonormal_vectors = gramschmidt(vectors)
%
%   INPUT ===========================================================
%
%   vectors (numeric matrix)
%   A matrix where each column represents a vector in R^n space.
%   Example: [1 2; 2 4]
%
%   OUTPUT ==========================================================
%
%   orthonormal_vectors (numeric matrix)
%   A matrix where each column is an orthonormal vector resulting from 
%   the Gram-Schmidt orthogonalization of the input vectors.
%
%   PROCESS =========================================================
%
%   The function works by iterating through each vector and subtracting 
%   its projection onto all previous orthonormal vectors, and then 
%   normalizing the result.
%
%   AUTHOR ==========================================================
%
%   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
%
%   =================================================================


% Get the number of vectors
num_vectors = size(vectors, 2);

% Initialize the matrix to store the orthonormal vectors
orthonormal_vectors = zeros(size(vectors));

% Apply the Gram-Schmidt orthogonalization process
for i = 1:num_vectors
    % Start with the current vector
    new_vector = vectors(:, i);

    % Subtract the projection onto the previous vectors
    for j = 1:i-1
        projection = dot(vectors(:, i), orthonormal_vectors(:, j)) / norm(orthonormal_vectors(:, j))^2;
        new_vector = new_vector - projection * orthonormal_vectors(:, j);
    end

    % Normalize the new vector
    orthonormal_vectors(:, i) = new_vector / norm(new_vector);
end
end
   