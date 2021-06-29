function m_dBm_grad = Watt_to_dBm_grad( m_Watt ) 
% The m_dBm(i,j) is the derivative of the scalar function
% f(x)=Watt_to_dBm(x) with respect to x and evaluated at m_Watt(i,j). 

m_dBm_grad = 10./m_Watt./log(10);

end