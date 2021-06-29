function m_dBm = Watt_to_dBm( m_Watt ) 

m_mWatt = m_Watt*1e3;
m_dBm = 10*log10( m_mWatt );

end