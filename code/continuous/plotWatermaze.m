function plotWatermaze( fileName, posXPlatform, posYPlatform )

%load
fid = fopen(fileName,'r');

allRecords=[];
recId=1;
while feof(fid) ~= 1
    line = fgetl(fid);
    
    rem=line;
    newRec=[];
    for i=1:9
        [token,rem] = strtok(rem, ',');
        if i==3
            %is in HEX
            newRec(i) = hex2dec( token );
        else
            newRec(i) = str2num( token );
        end
    end
    allRecords(recId, : ) = newRec;
    recId=recId+1;
end
fclose(fid);

%visualize
trials = unique(allRecords(:,1));

for i=2:length(trials)
    subplot(5,5,i-1);
    trial = trials(i);
    
    trialRecords =  allRecords( find(allRecords(:,1)==trial),:);
   
    %6, 8 is pos
    %7 is hight 
    
    plot( trialRecords(:,6), trialRecords(:,8),'x-','LineWidth',1); 
    hold on
    
    if trials(i)<21
        ylabel(['trial ' num2str(trials(i))]);
        plot( posXPlatform, posYPlatform, 'd', 'color','r', 'LineWidth',2);
    else
        ylabel(['probe trial ']);
    end
    xlim([-7 7]);
    ylim([-7 7]);
    hold off
end