function f_val = f_exp(q,p,coef)
% coef =sigma^2
f_val =1/(sqrt(2*pi*coef))* exp(-1/(2*coef)*norm(q -p)^2);