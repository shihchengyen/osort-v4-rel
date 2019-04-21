%
% plots details about a single cluster (waveforms, ISIs, powerspectra)
%
function plotSingleCluster(allSpikes, allTimestamps, assigned, label, clusterNr, spikeColor,timestampInclude)

if nargin==6
   timestampInclude=[];
end

spikeLength=size(allSpikes,2);

cluNr=clusterNr;
spikesToDraw = allSpikes( find(assigned==cluNr),:);
timestamps = allTimestamps( find(assigned==cluNr) );

%for plotting/figures -- plot less waveforms than there are (to reduce file size)
%spikesToDraw=spikesToDraw(1:10:size(spikesToDraw,1),:);

[~,~,L_R]=computeLratio( allSpikes ,assigned,cluNr,2); % compute L-ratio

ISI_fitWindows = [200 700]; %fit gamma function and plot of ISI from 0 to this value.

[f,Pxxn,tvect,Cxx,edges1,n1,yGam1,edges2,n2,yGam2,mGlobal,m1,m2,percentageBelow,CV]  = getStatsForCluster(spikesToDraw, timestamps, ISI_fitWindows);

nrSpikes = size(spikesToDraw,1);
fEstimate= nrSpikes / ((timestamps(end)-timestamps(1))/1e6);

manualPlotMode = 0;
if manualPlotMode
    AdBitVolts= 0.000000030518510385491027;  % set manual, just for plotting
    spikesToDraw = spikesToDraw * AdBitVolts * 1e6;  %now in uV
end

%% PDF
subplot(3,6,[5 6]);
spikePDFestimate(spikesToDraw(:,1:spikeLength));
set(gca,'XTickLabel',{});
title(['probability density function isolation distance=',num2str(L_R)])

%% waveforms
subplot(3,6,[1 2]);
% plot(1:spikeLength, spikesToDraw', spikeColor);
plot(1:spikeLength, spikesToDraw', 'color', spikeColor); % HM edit
xlim( [1 spikeLength] );
ylabel([label 'C' num2str(cluNr) ' n=' num2str(nrSpikes)]);
title(['Raw waveforms f=' num2str(fEstimate,3) 'Hz']);
hold on

avColor='b';

%Make sure the average does have a different color from the spikes
if spikeColor==avColor
    avColor='k';
end

plot(1:spikeLength, mean(spikesToDraw),avColor, 'linewidth', 2);

hold off
%ylim( [-1000 2000] );
set(gca,'XTickLabel',{});
%set(gca,'YTickLabel',{});

%xlim([50 250]);   % for plotting figures in papers
%ylim([-40 80]);

%power spectrum
subplot(3,6,[3 4]);
        
%convert to spiketrain
plot(f,Pxxn,'r','linewidth',2);        
xlim( [0 80] );
ylim( [0 max(Pxxn)*1.5+1] );  %+1 to prevent 0
xlabel('Hz');
ylabel('(spk/s)^2/Hz');

[isOk2]= checkPowerspectrum(Pxxn,f, 20.0, 100.0);  %check for peaks in powerspectrum in 20.0 ... 100.0 range
stat='yes';
if isOk2==false
    stat='no';
end

title(['Powerspectrum good=' stat]);

%std of waveform
subplot(3,6,[7 8])
S=std(spikesToDraw);
plot(1:256, S,'r','LineWidth',2');
line([95 95],[-3000 3000],'color','m');
xlim([1 256]);
ylim([min(S)*0.5 max(S)*1.5]);
title(['\sigma. \sigma(\sigma)=' num2str(std(S))]);
ylabel('\sigma');

%autocorrelation
subplot(3,6,9)
plot(tvect,Cxx,'r','LineWidth',2);
title('Autocorrelation (<10ms)');
ylabel('(spk/s)^2/Hz');
xlabel('[ms]');
xlim([1 10]);

%autocorrelation
subplot(3,6,10)
plot(tvect,Cxx,'r','LineWidth',2);
title('Autocorrelation (>10ms)');
ylabel('(spk/s)^2/Hz');
xlabel('[ms]');
xlim([10 80]);

%histograms
subplot(3,6,[13 14]);
bar(edges1,n1,'histc');


h=get(gca);
changeBARcolor(h.Children,[0 0 1]);

title(sprintf('ISI bin=1ms, mean=%.1f, below 3ms=%.2f%% [%.2f%%]',m1,percentageBelow(1),percentageBelow(3)));
%set(gca,'XTickLabel',{});
%set(gca,'YTickLabel',{});
ylabel( ['CV=' num2str(CV,3)] );
hold on
plot(edges1, yGam1, 'r','linewidth',2);
hold off
xlim( [0 ISI_fitWindows(1)] );
xlabel('ms');

subplot(3,6,[15 16]);
bar(edges2,n2,'histc');
title(sprintf('ISI bin=5ms, mean=%.1f, below 3ms=%.2f%%',m2,percentageBelow(2)));

h=get(gca);
changeBARcolor(h.Children,[0 0 1]);
% for i=1:10
% try
% set(h.Children(i),'FaceColor',[0 0 1]);
% set(h.Children(i),'EdgeColor',[0 0 1]);
% end
% end


hold on
plot(edges2, yGam2, 'r','linewidth',2);
hold off
xlim( [0 ISI_fitWindows(2)] );
xlabel('ms');

%% pdf, hist and timecourse kaminskij
subplot(3,6,[11 12]);
hist(spikesToDraw(:,95),50)
h=get(gca);
changeBARcolor(h.Children,[0 0 1]);
title('Waveform amp at alignment')
xlabel('amp [uV]');
ylabel('nr of spikes');

%% plot of firing rate and spike amplitude across time
%only taking into accout time form first and last detected spike ?

subplot(3,6,[17 18]);
l=(allTimestamps(end)-allTimestamps(1))/100;

for i=1:100
    % fire rate
  H(i)= length(find(timestamps>allTimestamps(1)+1+l*(i-1) ...
      & timestamps<allTimestamps(1)+l*(i)));
   % amp
  H1(i)= mean(spikesToDraw(find(timestamps>allTimestamps(1)+1+l*(i-1) ...
      & timestamps<allTimestamps(1)+l*(i)),95));
  
end

H=H/(l/1000000);
 
% for timestampInclude
if ~isempty( timestampInclude)
    
    timestampInclude2=(timestampInclude-allTimestamps(1))/1e6;
    
    for i=1:size(timestampInclude,1)
   
temp=round(timestampInclude2(i,1):timestampInclude2(i,2));

    %plot(temp,mean(H),'.g')
    
    area(temp,(max(H1)+1)*ones(length(temp),1),min(H1)-1,'FaceColor',[0.8 0.8 0.8],...
    'EdgeColor',[0.8 0.8 0.8]);
    hold on
    end
end

% ploting
[ax h1 h2]=plotyy(linspace(allTimestamps(1),allTimestamps(end),100)/1000000-allTimestamps(1)/1e6,H1,...
    linspace(allTimestamps(1),allTimestamps(end),100)/1000000-allTimestamps(1)/1e6,H);
try
ylim([min(H1) max(H1)]);
catch
 ylim([min(H1)-1 max(H1)]);   
end

set(h1,'LineW',2,'LineStyle','--');
set(h2,'LineW',2);

set(get(ax(1),'Ylabel'),'String','spike amplitude (\muV)') ;
set(get(ax(2),'Ylabel'),'String','firing rate (Hz)') ;

xlabel('time (sec)');
title('spike rate & amplitude during analyzed period');

axes(ax(1));
axis tight
axes(ax(2));
axis tight
 