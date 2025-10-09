% function [a, b] = mat2vecs(C) %deprecated
% %DECONSTRUCT_COMBINED_MATRIX Retrieves the two original vectors a and b 
% % from the combined matrix C.
% %
% %   [a, b] = DECONSTRUCT_COMBINED_MATRIX(C)
% %
% % Inputs:
% %   C - A matrix of size (n*m) x 2, created by combining vectors a and b.
% %
% % Outputs:
% %   a - The original numeric vector of size n.
% %   b - The original numeric vector of size m.
% %
% % IMPORTANT: This function assumes C was created using the 
% % 'combine_arrays_to_matrix' function provided previously, which implies
% % the rows vary according to 'a' first, then 'b'.
% 
% % 1. Determine the overall size of the matrix
% total_rows = size(C, 1);
% 
% % 2. Get the number of unique values in the second column (C(:, 2)).
% %    This corresponds to the length of the vector 'b', since 'b' changes
% %    less frequently in the 'C' matrix.
% %    The number of unique elements in C(:,2) is the length of 'b'.
% len_m = length(unique(C(:, 2))); % Length of vector 'b' (m)
% 
% % 3. Calculate the length of vector 'a' (n)
% %    Since total_rows = n * m, we have n = total_rows / m.
% len_n = total_rows / len_m;
% 
% % --- Reconstruct Vector a ---
% % Vector 'a' is in the first column (C(:, 1)). It repeats 'm' times.
% % We only need the first 'n' unique, consecutive values.
% % We can reshape the column into an (n x m) matrix and take the first column, 
% % or just take the first 'n' elements, since they should represent the vector 'a' 
% % before it starts repeating.
% A_col = C(:, 1);
% a = A_col(1:len_n)'; % Take the first 'n' elements and ensure it's a row vector
% 
% % --- Reconstruct Vector b ---
% % Vector 'b' is in the second column (C(:, 2)).
% % This column repeats the entire vector 'a' * len_m times.
% % We can reshape this column into an (n x m) matrix.
% B_col = C(:, 2);
% B_matrix = reshape(B_col, len_n, len_m);
% 
% % Now, 'b' is represented by the *first row* or *last row* of B_matrix.
% % For simplicity, since all rows should be identical to 'b', we extract the first row.
% b = B_matrix(1, :); % Extract the first row as vector 'b'
% 
% % If you want them as column vectors instead of row vectors (as they often are in MATLAB):
% % a = a';
% % b = b';
% 
% end