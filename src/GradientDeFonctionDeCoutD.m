function result = GradientDeFonctionDeCoutD(t, y, x)
    d = x(1);
    a = x(2);
    s = x(3);
    a2 = a.^2;
    s2 = ((1+s^2)^2)*2;
    s3 = s2./2;
    result = sum(4*a2/s2.*(-y.*(t-d).*exp(-(t-d).^2./s2)+a2.*(t-d).*exp(-(t-d).^2./s3)));
end