function result = GradientDeFonctionDeCout(t, sig_noisy, xk)
    result = [GradientDeFonctionDeCoutD(t,sig_noisy,xk); GradientDeFonctionDeCoutA(t,sig_noisy,xk); GradientDeFonctionDeCoutS(t,sig_noisy,xk)];
end