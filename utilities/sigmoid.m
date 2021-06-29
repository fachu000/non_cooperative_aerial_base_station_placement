function m_out = sigmoid(m_in)

% Entrywise application of the sigmoid function  @(x) exp(x)./(exp(x)+1)
% with enhancements for numerical stability

m_out =1./(exp(-m_in)+1);

end