function [x] = w2dbm(power)
    digits(32);
    x = vpa(10 * log10(power*1000));
end