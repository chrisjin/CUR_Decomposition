function [result, c, u ,r] = cur_solver(mat, columnarr, rowarr)
cn = length(columnarr);
cr = length(rowarr);
width = cn;

matsquare = mat.^2;

tnorm = sum(sum(matsquare));
cnorm = sum(matsquare, 1);
rnorm = sum(matsquare, 2);

columns = mat(:, columnarr);
qcol = cnorm(columnarr)/tnorm;
for i = 1:width
    columns(:, i) = columns(:, i) / sqrt(width * qcol(i));
end

rows = mat(rowarr, :);
qrow = rnorm(rowarr)/tnorm;
for i = 1:width
    rows(i, :) = rows(i, :) / sqrt(width * qrow(i));
end

c = columns;
r = rows;

w = zeros(width, width);
for i = 1:width
    for j = 1:width
        w(i, j) = mat(rowarr(i), columnarr(j));
    end
end

[x, lam, y] = svd(w);
lamplus = pinv(lam);
u = y*lamplus*lamplus*x';

result = c*u*r;





