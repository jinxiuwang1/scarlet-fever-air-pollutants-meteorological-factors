function [yy, dy] = B_spline_der(x, y, h)
% Sort the data (ordered data required for B-spline)
[x_sorted, sort_idx] = sort(x);
y_sorted = y(sort_idx);

% Set the B-spline order (k=4 for cubic spline)
k = 4;

% Generate the B-spline interpolation function
spline = spapi(k, x_sorted, y_sorted);

% Compute the first derivative (key step)
spline_der = fnder(spline, 1);  % 1 indicates first derivative

% Generate dense points for plotting the curve
xx = min(x_sorted):h:max(x_sorted);

% Compute the original B-spline interpolation values
yy = fnval(spline, xx);

% Compute the first derivative values
dy = fnval(spline_der, xx);
end