function [color_img] = plot_all_contours( img, array_of_ellipses, plot_bool)
%%
% Author: Camilo Aguilar
% plots all the contour on top of the image img
% if plot_bool=1, the function plots an image in Matlab's GUI too
%%
    img = uint8(img);
    if(size(img,3)==1)
        color_img = cat(3, img, img, img);
    else
        color_img=img;
    end
    
    for ellipse=1:length(array_of_ellipses)
               
        %POINT_NUMBER = array_of_ellipses(ellipse).p_n;  
            cx = round(array_of_ellipses(ellipse).rx);
            cy = round(array_of_ellipses(ellipse).ry);
            POINT_NUMBER = length(cx);
            r_value = max(max(color_img(:,:,1)));
            if(r_value == 0)
                r_value = 1;
            end
            for k=1:POINT_NUMBER-1
                color_img(cy(k),cx(k),1) = r_value;
                color_img(:,:,1) = conect_line( color_img(:,:,1),cx(k),cy(k),cx(k+1),cy(k+1),r_value);
            end
            color_img(:,:,1) = conect_line( color_img(:,:,1),cx(POINT_NUMBER),cy(POINT_NUMBER),cx(1),cy(1),r_value);
        
              
    end
    if(exist('plot_bool','var'))
        imagesc(color_img);
        hold on;
         for ellipse=1:length(array_of_ellipses)
            x = mean(array_of_ellipses(ellipse).rx);
            y = mean(array_of_ellipses(ellipse).ry);
            text(x,y,num2str(ellipse),'fontsize',10) 
            try
                text(x+5,y+5,num2str(array_of_ellipses(ellipse).E),'fontsize',10) 
            catch 
                %do nothing
            end

         end
        hold off;
         title('Detected Objects');
    end
    
    
    

   
end

