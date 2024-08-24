function [rank] = count_rank(vector, factor)
    % This function counts the number of eigenvalues distinctly larger than zero

    % Initialize the count of distinct eigenvalues
    if vector(1) == 0
        error('Invalid vector')
    else
        rank = 1;
    end

    % Loop through the vector to count the eigenvalues distinctly larger than zero
    for i = 1:length(vector) - 1
        if vector(i + 1) / vector(i) > factor
            rank = rank + 1;
        else
            % Once we find an element that doesn't meet the condition, we break out of the loop
            break;
        end
    end

end