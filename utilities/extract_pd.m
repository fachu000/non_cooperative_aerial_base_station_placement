

function pd=extract_pd(pfa,d0,d1)
% PD = EXTRACT_PD(PFA,D0,D1)
% pd is the probability of detection computed row-wise, i.e, every row is
% treated independently. pd is a column vector with the same number of rows
% as d0.

% compute detection
pd=zeros(size(d0,1),1);
for d=1:size(d0,1)
	thr=prctile(d0(d,:),100*(1-pfa));
	
	%pfa=sum(d0(d,:)>thr)/size(d1,2)
	
	pd(d)=sum(d1(d,:)>thr)/size(d1,2);
	
end

return

