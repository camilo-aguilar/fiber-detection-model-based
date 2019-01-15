function [list_of_contours] = get_list_of_contours(c)
%%
% Author: Camilo Aguilar
% Gets a level sets and returns a list of contour objects
%%

idx = 1;
count = 1;
if nargout > 1, contour_pt = []; end
dim = size(c,1); 
list_of_contours(1).rx = [];
list_of_contours(1).ry = [];

while idx < size(c,2)
    n = c(2,idx);
    if(n > 30)
        list_of_contours(count).rx = c(1,idx+1:idx+n);
        list_of_contours(count).ry = c(2,idx+1:idx+n);    
        count = count+1;
    end
    
    idx = idx+n+1;
end