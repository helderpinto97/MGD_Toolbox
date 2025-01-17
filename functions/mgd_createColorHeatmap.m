function mgd_createColorHeatmap(cellArray)
    % Step 1: Find the most frequent vectors
    mostFrequentVectors = findMostFrequentVectors(cellArray);

    % Step 2: Initialize a heatmap matrix
    heatmapMatrix = zeros(size(cellArray));

    % Step 3: Compare each cell entry with the most frequent vector for the column
    for col = 1:size(cellArray, 2)
        % Get the most frequent vector for this column
        frequentVector = sort(mostFrequentVectors{col}); % Sort for comparison
        
        % Compare each cell in the column with the most frequent vector
        for row = 1:size(cellArray, 1)
            currentVector = sort(cellArray{row, col}); % Sort for comparison
            if isequal(currentVector, frequentVector)
                heatmapMatrix(row, col) = 1; % Match
            else
                heatmapMatrix(row, col) = 0; % Mismatch
            end
        end
    end

    % Step 4: Display the heatmap
    figure;
    imagesc(heatmapMatrix); % Display heatmap
    colormap(jet); % Use a colorful colormap
    colorbar('Ticks', [0, 1], 'TickLabels', {'Mismatch', 'Match'}); % Add legend
    xlabel('Columns');
    ylabel('Rows');
    title('Color Heatmap of Matches with Most Frequent Vectors');
end

