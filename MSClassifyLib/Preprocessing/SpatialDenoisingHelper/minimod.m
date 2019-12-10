function fluxLimiter = minimod(a,b)
% Used e.g. in the shock-Filter, if a and b are two approximations of a
% simple finite difference minimod chooses the absolute-wise smaller one or
% 0 in case of different signs

fluxLimiter = sign(a).*min(abs(a),abs(b));
fluxLimiter(a.*b<=0) = 0;

end