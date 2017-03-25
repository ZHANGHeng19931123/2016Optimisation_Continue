function result = GradientDeFonctionDeCoutS(t, y, x)
    d = x(1);
    a = x(2);
    s = x(3);
    a2 = a^2;
    s2 = ((1+s^2)^2)*2;
    result = -4.*sum(a2.*s.*(y-a2.*exp(-(t-d).^2./s2)).*exp(-(t-d).^2./s2).*((t-d).^2)./((1+s.^2).^3));
end