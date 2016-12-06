function [solution]=cwm(p,p2,rs, thresnum, threscon)
timerVal = tic;
%% Parameters
[I,J] = size(p); % I1 is number of forecasters, J1 is number of events
% Score coeffiecients:
a = 101;
b = 100/log(a);

%% calculate number of forecasting
cnt = zeros(I,1); % count of the number of events forecasted by each forecaster
num = zeros(J,1); % number of forecasts received by each event
ind = zeros(I,J); % indicator if forecaster i forecasted event j
for i=2:I
    for j=2:J
        if isnan(p(i,j))==0
            cnt(i,1)=cnt(i,1)+1;
            ind(i,j)=1;
            num(j,1)=num(j,1)+1;
        end
    end
end
%% calculate m
m = zeros(j,1);
for j=2:J
    temp = 0;
    for i=2:I
        if isnan(p(i,j)) == 0
            temp = temp + p(i,j);
        end
    end
    m(j,1) = temp/sum(ind(:,j));
end

%% calculate score of jth event with ith forecaster
SJ = zeros(J,1); % score for jth event
for j=2:J
    SJ(j,1)=b*log(a-rs(j,1)*m(j,1)-(1-rs(j,1))*(100-m(j,1)));
end

%% calculate mi
mi = zeros(i,j);
for i=2:I
	for j=2:J
        temp = 0;
        ct = 0;
        for ii=2:I
            if (isnan(p(i,j)) == 0) && (ii ~= i)
                temp = temp + p(i,j);
                ct = ct + 1;
            end
        end
        mi(i,j) = temp/ct;
	end
end

%% calculate score of jth event without ith forecaster
SIJ = zeros(I,J); % score for ith event without jth forecaster
for j=2:J
    for i=2:I
        SIJ(i,j)=b*log(a-rs(j,1)*mi(i,j)-(1-rs(j,1))*(100-mi(i,j)));
    end
end

%% calculate their contribution
contribution = zeros(I,1); % contribution of ith forecaster
for i=2:I
    if cnt(i,1) >= thresnum
        temp = 0;
        ct = 0;
        for j=2:J
            if  isnan(p(i,j)) == 0
                temp = temp + SJ(j,1) - SIJ(i,j);
                ct = ct + 1;
            end
        end
    end
    contribution(i,1) = temp / ct;
end

%% crowds forecasting assigned probability
[I2,J2] = size(p2);
prob = zeros(J2,1);
for j=2:J2
    temp = 0;
	tempcontri = 0;
    for i=2:I2
        if (isnan(p2(i,j)) == 0) && (contribution(i,1)>=threscon)
            temp = temp + p2(i,j) * contribution(i,1);
            tempcontri = tempcontri + contribution(i,1);
        end
    end
	prob(j,1) = temp / tempcontri;
end

%% Individual scores:
s = zeros(I2,J2);
S = zeros(I2,1);
cnt2 = zeros(I2,1);
for i=2:I2
    for j=2:J2
        if isnan(p2(i,j)) == 0
            s(i,j)=b*log(a-rs(j,1)*p2(i,j)-(1-rs(j,1))*(100-p2(i,j)));
            cnt2(i,1) = cnt2(i,1) + 1;
        end
    end
    S(i,1) = sum(s(i,:))/cnt2(i,1); 
end

%% Calculate the aggregating scores
sc = zeros(J2,1);
for j=2:J2
    sc(j,1)=b*log(a-rs(j,1)*prob(j,1)-(1-rs(j,1))*(100-prob(j,1)));
end
score = mean(sc);
solution.score = score;
solution.cnt = cnt;
solution.num = num;
solution.m = m;
solution.mi = mi;
solution.SJ = SJ;
solution.SIJ = SIJ;
solution.contribution = contribution;
solution.prob = prob;
solution.sc = sc;
solution.t = toc(timerVal);
solution.S = S;
solution.cnt2 = cnt2;

