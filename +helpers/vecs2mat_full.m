function C = vecs2mat_full(a,b)            

% exhaustive combination
[B_grid, A_grid] = ndgrid(b, a);

% 2. Reshape the grids into a single column vector
A_col = A_grid(:);
B_col = B_grid(:);

% 3. Concatenate the column vectors side-by-side to form the final matrix
C = [A_col, B_col];
end