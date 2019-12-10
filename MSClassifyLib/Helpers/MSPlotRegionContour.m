function[contourPlot] = MSPlotRegionContour(region, lineWidth, lineColor,ax)
% This function plots the contour of a specified region to the current
% figure. The contour line is smoothed by convolution with a gaussian
% kernel. 
%
% Parameters:
% region     - Region whose contour will be plotted as a matrix
% lineWidth  - Line width of the contour plot. If none is given, a default
%              value of 2 is chosen.
% lineColor  - Color of the contour plot. The default color is red.
% interpolationFactor - Value for resizing the region (is needed for the
% m/z visualizing GUI)
    if(exist('region','var'))
        % Convert region, if it is not of type single or double
        if( ~any(ismember(class(region),{'single', 'double'})))
            region = single(logical(region));
        end
        
        %Smoothen the contour line
        filter = fspecial('gaussian',2,1);
        region = conv2(region,filter); 
 
        minimum = min(region(region > 0));
        if(exist('ax','var'))
            [~, contourPlot] = contour(ax,region,[minimum,minimum]);
        else
            [~, contourPlot] = contour(region,[minimum,minimum]);
        end

        %Set contour plot properties
        if(exist('lineWidth','var') )
            set(contourPlot,'LineWidth',lineWidth) 
        else
            set(contourPlot,'LineWidth',2); 
        end
        
        if(exist('lineColor','var') )
            set(contourPlot,'LineColor',lineColor);
        else
            set(contourPlot,'LineColor','red');
        end

    else
        error('No region was specified.');
    end
end
