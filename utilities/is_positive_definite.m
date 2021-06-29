function d = is_positive_definite( C )

[~,p] = chol(C);
d = (p == 0);

end