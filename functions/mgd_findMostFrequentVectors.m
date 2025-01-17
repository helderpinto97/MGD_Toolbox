function mostFrequentVectors = mgd_findMostFrequentVectors(cellArray)
    % Initialize output array
    mostFrequentVectors = cell(1, size(cellArray, 2));
    
    % Iterate over each column of the cell array
    for col = 1:size(cellArray, 2)
        % Get the column vectors
        columnData = cellArray(:, col);
        
        % Sort each vector to remove the effect of the order
        sortedVectors = cellfun(@(x) sort(x), columnData, 'UniformOutput', false);
        
        % Convert sorted vectors to strings for comparison
        stringVectors = cellfun(@(x) mat2str(x), sortedVectors, 'UniformOutput', false);
        
        % Find unique string vectors and their frequencies
        [uniqueVectors, ~, idx] = unique(stringVectors, 'stable');
        
        % Count the frequency of each unique vector
        frequency = histcounts(idx, 1:length(uniqueVectors));
        
        % Find the index of the most frequent vector
        [~, maxIndex] = max(frequency);
        
        % Convert the most frequent string vector back to a numeric vector
        mostFrequentVectors{col} = str2num(uniqueVectors{maxIndex}); %#ok<ST2NM>
    end
end
