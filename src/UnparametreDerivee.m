function result = UnparametreDerivee(t, y, d)
    s = 2;
    a = 2;
    a2 = a^2;
    s2 = ((1+s^2)^2)*2;
    s3 = s2./2;
    result = zeros(size(d));
    for i = 1:size(d,2)
        result(i) = sum(4*a2/s2.*(-y.*(t-d(i)).*exp(-(t-d(i)).^2./s2)+a2.*(t-d(i)).*exp(-(t-d(i)).^2./s3)));
    end
end