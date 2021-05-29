function [mcf,mcfll,mcful,time]=mcfun(events,censor)

censor=horzcat(censor{:});

events=horzcat(events{:});
events(isnan(events))=[];

time=unique(sort(events));
atrisk=bsxfun(@ge,censor',time);

eventsPerTime=sum(bsxfun(@eq,sort(events),time'),2);
mcf=cumsum(eventsPerTime'./sum(atrisk));

totalatrisk=sum(atrisk);
mcfvar=cumsum(1./totalatrisk.^2.*((1-1./totalatrisk).^2+(totalatrisk-1).*1./totalatrisk.^2));
mcfll=mcf./exp(norminv(0.975).*sqrt(mcfvar)./mcf); %https://doi.org/10.1198/004017008000000091
mcful=mcf.*exp(norminv(0.975).*sqrt(mcfvar)./mcf);

mcf=reshape(repmat(mcf,2,1),1,[]);
mcf=mcf(1:end-1);
mcfll=reshape(repmat(mcfll,2,1),1,[]);
mcfll=mcfll(1:end-1);
mcful=reshape(repmat(mcful,2,1),1,[]);
mcful=mcful(1:end-1);

mcf(end+1)=mcf(end);
mcfll(end+1)=mcfll(end);
mcful(end+1)=mcful(end);

time=reshape(repmat(time,2,1),1,[]);
time=time(2:end);
time(end+1)=max(censor);