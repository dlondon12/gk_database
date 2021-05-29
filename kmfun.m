function [km,kmll,kmul,time]=kmfun(events,censor)

events=cellfun(@(x) x(1),events);
dbfkp=events;
events(isnan(events))=[];
time=unique(sort(events));
censormat=horzcat(censor{:});


eventsPerTime=sum(bsxfun(@eq,sort(events),time'),2);
cens=min([dbfkp',censormat'],[],2,'omitnan');
atrisk=bsxfun(@ge,cens,time);
totalatrisk=sum(atrisk);
km=cumprod(1-eventsPerTime'./totalatrisk);

kmvar=1./(log(km).^2).*cumsum(eventsPerTime'./(totalatrisk.*(totalatrisk-eventsPerTime')));
cplus=log(-log(km))+norminv(0.975)*sqrt(kmvar); %exponential greenwood confidence interval, Kalbfleisch and Prentice 1980
cminus=log(-log(km))-norminv(0.975)*sqrt(kmvar);
kmul=exp(-exp(cminus));
kmll=exp(-exp(cplus));
km=[1,km];
kmll=[1,kmll];
kmul=[1,kmul];

km=reshape(repmat(km,2,1),1,[]);
km=km(1:end-1);
kmll=reshape(repmat(kmll,2,1),1,[]);
kmll=kmll(1:end-1);
kmul=reshape(repmat(kmul,2,1),1,[]);
kmul=kmul(1:end-1);

km(end+1)=km(end);
kmll(end+1)=kmll(end);
kmul(end+1)=kmul(end);

time=reshape(repmat([0,time],2,1),1,[]);
time=time(2:end);
time(end+1)=max(censormat);
end