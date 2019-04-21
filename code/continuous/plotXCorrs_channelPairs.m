%
% for pairs of electrodes, plot
%
%inds is, one per row, from/to/panelNr
%plotMode=1 is coeff, =2 is biased(raw)
%
%urut/MPI/sept11
function plotXCorrs_channelPairs(dataSpikesChannel,channels, inds, Fs, maxlag, plotRange,sessionID,blockNrUse, CARcorrection, plotMode, plotIteration)
if nargin<9
    CARcorrection=[];
end
if nargin<10
    plotMode=2;
end
if nargin<11
    plotIteration=1;
end

%%
for k=1:size(inds,1)

    ax(k)=subplot( 4,4,inds(k,3));

    if plotIteration>1
        hold on;
        col='r';
    else
        col='b';
    end
    
    disp(['from /to' num2str(inds(k,1:2))]);

    d1=dataSpikesChannel{inds(k,1)}(plotRange);
    d2=dataSpikesChannel{inds(k,2)}(plotRange);
    
    if ~isempty(CARcorrection)
        d1 = d1-CARcorrection;
        d2 = d2-CARcorrection;
        pstrCAR='CARcorr=ON';
    else
        pstrCAR='CARcorr=OFF';
    end
    
    if plotMode==1
        c = xcorr( d1, d2, maxlag,'coeff');
    else
        c = xcorr( d1, d2, maxlag,'biased');
    end
    
    t=[-maxlag:maxlag]/Fs*1000;


    plot(t,c, col, 'linewidth',2)
 %   xlim([-4 4]);
    
    if k==1
        pstr=[sessionID '-b' num2str(blockNrUse)];
    else
        pstr='';
    end
    
    title([pstrCAR ' ' pstr ' xcorr Ch' num2str(channels(inds(k,1))) '/' num2str(channels(inds(k,2)))],'interpreter','none');
    xlabel(['ms']);
    ylabel(['max=' num2str(max(c)) ]);
    if plotMode==1
    ylim([-0.5 1.0]);
    end
    hold off
end
linkaxes(ax,'xy');
xlim([-4 4]);
