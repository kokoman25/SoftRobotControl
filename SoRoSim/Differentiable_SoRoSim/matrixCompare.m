%compares 2 matricies of same size
function matrixCompare(A,B)
    % Compute absolute difference matrix
    diff_matrix = abs(A - B);
    
    % Compute relative difference (ignoring values < 1e-6)
    valid_mask = max(abs(A), abs(B)) > 1e-6;  % Ignore small values
    relative_diff = zeros(size(A));  % Initialize with zeros
    relative_diff(valid_mask) = diff_matrix(valid_mask) ./ max(abs(A(valid_mask)), abs(B(valid_mask)));
    
    % Find the maximum relative difference and its index
    [max_rel_value, linear_idx] = max(relative_diff(:));
    
    % Convert linear index to row and column indices
    [row, col] = ind2sub(size(A), linear_idx);
    
    % Display the results
    fprintf('Max absolute difference: %e at row %d, column %d\n', diff_matrix(row, col), row, col);
    fprintf('Max relative difference: %.2f%% at row %d, column %d\n', max_rel_value * 100, row, col);
end