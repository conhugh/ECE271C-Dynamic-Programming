function features = extract_features(board)

    heights = boardHeight(board);
%     max_height = max(heights);
%     d2top = 15 - max_height;
%     height_diffs = heights(2:length(heights))-heights(1:length(heights)-1);
%     avg_diff = mean(abs(height_diffs));
%     avg_height = mean(heights);
%     avg_h2top = 15 - avg_height;
%     stdev = std(heights);
    
    square_heights = sum(heights.^2)/50;
    
    holes = sum(heights - sum(board,1));

    %features = [d2top; avg_diff; avg_h2top; stdev; holes];
    features = [square_heights; holes];
    
end

