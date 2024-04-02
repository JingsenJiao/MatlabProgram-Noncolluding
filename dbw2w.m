function [power] = dbw2w(x)
    digits(32);
    power = vpa(10^(x/10));
end