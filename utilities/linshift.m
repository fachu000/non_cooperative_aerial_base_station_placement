function Mout = linshift(Min,r)

if length(r)==1
	Mout= linshift_onedim(Min,r);
elseif length(r)==2
	r1 = r(1);
	Mout = linshift_onedim(Min,r1);
	r2 = r(2);
	Mout = linshift_onedim(Mout',r2)';
	
else
	error('size of R not implemented');
end

end


function Mout = linshift_onedim(Min,r)
	Mout = circshift(Min,r);
	if r>=0
		Mout(1:r,:)=0;
	else
		Mout(r+end+1:end,:)=0;
	end
end