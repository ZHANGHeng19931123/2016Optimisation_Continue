function result = UnparametreDeriveeSecond(t, y, d)
    s = 2;
    a = 2;
    a2 = a^2;
    s2 = ((1+s^2)^2)*2;
    s3 = (1+s^2)^2;
    result = zeros(size(d));
    for i = 1:size(d,2)
        p1 = -(t-d(i)).^2./s2;
        p2 = -(t-d(i)).^2./s3;
        result(i) = sum(2*a2/s3.*(((t-d(i)).^2./s3).*(-y.*exp(p1)+2.*a2.*exp(p2))+y.*exp(p1)-a2.*exp(p2)));
    end
end