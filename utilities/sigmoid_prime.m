function m_out = sigmoid_prime(m_in)

% Entrywise application of the first-order derivative of the sigmoid 
% function  sigmoid = @(x) exp(x)./(exp(x)+1) with enhancements for
% numerical stability 

m_out = sigmoid(m_in).*(1-sigmoid(m_in));

end