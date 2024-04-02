function [x] = w2dbw(power)
    digits(32);
    x = vpa(10 * log10(power));
end