function yy=B_spline(x,y,h)

% Sort the data (ordered data required for B-spline)
[x_sorted, sort_idx] = sort(x);
y_sorted = y(sort_idx);

% Set the B-spline order (k=4 for cubic spline, k=3 for quadratic)
k = 4;

% Generate the B-spline interpolation function
spline = spapi(k, x_sorted, y_sorted);

% Generate dense points for plotting the curve
xx=min(x_sorted):h:max(x_sorted);
yy = fnval(spline, xx);