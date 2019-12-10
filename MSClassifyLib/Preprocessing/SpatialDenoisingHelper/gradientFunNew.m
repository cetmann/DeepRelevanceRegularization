function [gx,gy] = gradientFunNew(im)
% Robust estimation of the gradient (used e.g. in the shock filter)

gx = minimod(diffDir(im,'x+'),diffDir(im,'x-'));
gy = minimod(diffDir(im,'y+'),diffDir(im,'y-'));

end

