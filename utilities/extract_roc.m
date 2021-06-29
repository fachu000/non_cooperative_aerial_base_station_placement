function  [pfa,pd]=extract_roc(d0,d1)
% [pfa,pd]=extract_roc(d0,d1)
%

niter=size(d0,2);
npoints=120;
pfa=zeros(size(d0,1),npoints);
pd=zeros(size(d0,1),npoints);
for d=1:size(d0,1)
	d01=sort([d0(d,:),d1(d,:)]);
	thresh= d01(round(linspace(1,length(d01),npoints)));
	for k=1:length(thresh)
		pfa(d,k)=sum(d0(d,:)>thresh(k))/niter;
		pd(d,k)=sum(d1(d,:)>thresh(k))/niter;
	end
end

end