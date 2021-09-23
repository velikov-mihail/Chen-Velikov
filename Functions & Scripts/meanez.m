function mu = meanez(x,dim)

if nargin == 1
    dim = 1;
end

x = cropnan(x);

mu = mean(x,dim);