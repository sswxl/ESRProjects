function [datamat, ind1, ind2] = create_datamat(data, knots, spline_order)

    datamat = zeros(length(data), length(knots));
    ind1 = zeros(length(data), length(knots));
    ind2 = zeros(length(data), length(knots));
    
    for j = 1:length(knots)%这是一次样条，二次样条在下面
        datamat(:, j) = (data > knots(j)) .* (data - knots(j)).^spline_order;
        ind1(:, j) = (data > knots(j));
    end
    
    % 
    % datamat = zeros(length(data), length(knots) + spline_order);
    % ind1 = zeros(length(data), length(knots) + spline_order);
    % ind2 = zeros(length(data), length(knots) + spline_order);
    % for i = 1:spline_order
    %     datamat(:, i) = data.^i;
    %     ind1(:, i) = i * data.^(i - 1);
    %     ind2(:, i) = i * (i - 1) * data.^(i - 1);
    % end
    % for j = 1:length(knots)
    %     datamat(:, j + spline_order) = (data > knots(j)) .* (data - knots(j)).^spline_order;
    %     ind1(:, j + spline_order) = (data > knots(j));
    % end
    knots = linspace(0, 1, 100 + 1);
    xstar = linspace(0, 1, 101);