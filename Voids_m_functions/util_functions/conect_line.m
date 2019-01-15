function [ tmpmask ] = conect_line( tmpmask,xo,yo,xf,yf,new_c)
    if(xf == xo)
        while(yo ~= yf)
            step = (yf-yo)/abs(yf-yo);
			yo = yo + step;
			tmpmask(yo,xo) = new_c;
        end
        tmpmask(yf,xo) = new_c;
        return;
    end

	m = (yf-yo)/(xf-xo);
	b = yf - xf*m;
	
    while(xo ~= xf || yo ~= yf)
        tempX = xo + (xf-xo)/abs(xf-xo);
        tempY = floor(m*tempX + b + 0.5);
        tmpmask(tempY,tempX) = new_c;
		
        while(tempY - yo > 1)
            yo = yo + 1;
            tmpmask(yo,tempX) = new_c;
        end
		
		
        while(tempY - yo < -1)
			yo = yo - 1;
			tmpmask(yo,tempX) = new_c;
        end
		
		xo = tempX;
		yo = tempY;
    end

end

