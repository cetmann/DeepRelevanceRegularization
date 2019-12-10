function dq = diffDir(im,direction)
% Four basic difference schemes for images (used e.g. in the shock filter)

switch direction
  case 'x+'
    dq = diff([im(:,1) im],1,2);
  case 'x-'
    dq = -fliplr(diff(fliplr([im(:,1) im]),1,2));
  case 'y+'
    dq = diff([im(1,:); im],1,1);
  case 'y-'
    dq = -flipud(diff(flipud([im(1,:); im]),1,1));
end

end
