function [r,c] = fedaCore(rows_arr,cols_arr)
%Function to assist on FEDA's execution
r1 = rows_arr(1); r2 = rows_arr(length(rows_arr));
c1 = cols_arr(1); c2 = cols_arr(length(cols_arr));

Lr = r2-(r2-r1+1)/2;
Lc = c2-(c2-c1+1)/2;

r(1,:) = r1:Lr; r(2,:) = r(1,:); r(3,:) = Lr+1:r2; r(4,:) = r(3,:);
c(1,:) = c1:Lc; c(2,:) = Lc+1:c2; c(3,:) = c(1,:); c(4,:) = c(2,:);
end

