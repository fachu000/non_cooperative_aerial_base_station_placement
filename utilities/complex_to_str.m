function str=complex_to_str(patt,z)
% complex to string
%    patt  can be
%        '%d'
%        '%.2f'
%        '%.4f'
%        etc
if real(z)==0
	if imag(z)==0
		str=0;
	else
		if imag(z)==1
			str=sprintf('j');
		elseif imag(z)==-1
			str=sprintf('-j');
		else
			if imag(z)>0
				str=sprintf(['j' patt],imag(z));
			else
				str=sprintf(['-j' patt],abs(imag(z)));
			end
		end
	end
else
	if imag(z)==0
		str=sprintf(patt,real(z));
	else
		if imag(z)==1
			str=sprintf([patt '+j'],real(z));
		elseif imag(z)==-1
			str=sprintf([patt '-j'],real(z));
		else
			if imag(z)>0
				str=sprintf([patt '+j' patt],real(z),imag(z));
			else
				str=sprintf([patt '-j' patt],real(z),abs(imag(z)));
			end
		end
	end
end

return