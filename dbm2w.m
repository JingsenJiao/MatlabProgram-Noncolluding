function [power] = dbm2w(x)
%     digits(32);
%     power = vpa(10^(-3) * 10^(x/10));
    power = 10^(-3) * 10^(x/10);
end