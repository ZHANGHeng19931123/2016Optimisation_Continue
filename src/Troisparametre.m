function result = Troisparametre(t, y, x)
    d = x(1);
    a = x(2);
    s = x(3);
    a2 = a.^2;
    s2 = ((1+s.^2).^2).*2;
    result = sum((y-exp(-(t-d).^2./s2).*a2).^2);
end