%
%plots and exports  xcorr figures
% kaminskij 19 Dec 2014
% kaminskij 6 jan 2015 when scaling factor is nan (DBS data)
function produceXcorrFigures(handles, outputpath,outputFormat, thresholdMethod,XcorrlengthinMS,binsize)


%make dir to store the figures
if exist(outputpath)==0
    mkdir(outputpath);
end

colors = defineClusterColors;
rawWaveforms = handles.newSpikesNegative;

%order of clusters
%first clusters are the biggest clusters
clusters=handles.useNegative;
clusters = flipud( clusters );

%if there are more than we have color definitions,crop.
if length(clusters)>length(colors)
    clusters=clusters(1:length(colors));
end

nrClusters=length(clusters);




[outputEnding,outputOption] = determineFigExportFormat( outputFormat );

%--
%--- significance test between all clusters
pairs=[];
pairsColor=[];
c=0;

%only for the 7 biggest to avoid lots of useless plots
if nrClusters>7
    nrClusters=7;
end
%paperPosition=[0 0 16 12];
%%

%setfullscreen(1);
figure
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16*1.4 8*1.4]);
set(gcf,'visible','off');
%set(gcf,'PaperUnits','inches','PaperPosition',paperPosition)
timevector=(1:binsize:(XcorrlengthinMS*2+1))-XcorrlengthinMS-1;

for i=1:nrClusters
    for j=1:nrClusters
        c=c+1;
        pairs(c,1:2)=[clusters(i) clusters(j)];
        pairsColor(c,1:2)=[i j];
        
        if  (i)<(j)
            
            
            timestamps = handles.allSpikesTimestampsNegative( find(handles.assignedClusterNegative==clusters(i)) );
            n1 = convertToSpiketrain(timestamps,binsize);
            timestamps = handles.allSpikesTimestampsNegative( find(handles.assignedClusterNegative==clusters(j)) );
            n2 = convertToSpiketrain(timestamps,binsize);
            
            % making trains the same length for xcorr coeff normalization
            if length(n1)>length(n2)
                n2(length(n1))=0;
            elseif length(n1)<length(n2)
                n1(length(n2))=0;
            end
            Cxx=xcorr(n1,n2,XcorrlengthinMS/binsize,'coeff');
            
            
            %     for t=1:XcorrlengthinMS  X=[n1;n2]; plot(X(:,:)')
            %         Cxxdiff(t)=Cxx(XcorrlengthinMS+1-t)-Cxx(XcorrlengthinMS+1+t);
            %     end
            
            subplot('position',[0.06+(i-1)*0.93/nrClusters,0.06+(j-1)*0.93/nrClusters,0.7/nrClusters,0.65/nrClusters])
            area(timevector, Cxx)
            
            if j==nrClusters
                if i==1
                    title(['Cluster ',num2str(clusters(i)),' crosscorrelation'])
                else
                    title(['Cluster ',num2str(clusters(i))])
                end
            end
            if i==1
                
                ylabel(['Cluster ',num2str(clusters(j))])
                
            end
            h=get(gca);
            changeBARcolor(h.Children,[0 0 1]);
            %       subplot('position',[0.05+(j-1)*0.93/nrClusters,0.05+(i-1)*0.93/nrClusters,0.8/nrClusters,0.75/nrClusters])
            %     area( Cxxdiff)
            %     xlabel('time (ms)')
            %     h=get(gca);
            % changeBARcolor(h.Children,[0 0 1]);
        end
        if  (i)==(j)
            
            subplot('position',[0.06+(i-1)*0.93/nrClusters,0.06+(j-1)*0.93/nrClusters,0.7/nrClusters,0.65/nrClusters])
            if isnan(handles.scalingFactor)
                spikePDFestimate( handles.allSpikesNegative( find(handles.assignedClusterNegative==clusters(i)),: ));
                
            else
                spikePDFestimate( handles.allSpikesNegative( find(handles.assignedClusterNegative==clusters(i)),: )*handles.scalingFactor*1e6);
            end
            if i==1
                ylabel(['Cluster ',num2str(clusters(j))])
            end
            if j==nrClusters
                title(['Cluster ',num2str(clusters(i)),' Binsize=',num2str(binsize),' ms'])
            end
        end
    end
end
%set(gcf,'PaperPositionMode','auto')
%print(gcf, '-dtiff', [outputpath handles.prefix handles.from ' _Xcorr_THM_' num2str(thresholdMethod)]);
print(gcf, outputOption, [outputpath handles.prefix handles.from ' _Xcorr_THM_' num2str(thresholdMethod) outputEnding ]);
close(gcf);
%%



