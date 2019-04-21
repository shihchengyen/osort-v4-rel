function [REF REF2]=findREFelectrode(EEG)

EEG = pop_rmbase( EEG, [0      EEG.xmax*1000]);
for i=1:(size(EEG.data,1)-1)
    
   channelStd(i)=std(EEG.data(i,:)); 

end


figure;
a=pop_spectopo(EEG, 1, [0      834490.7237], 'EEG' , 'percent', 15, 'freqrange',[2 25],'electrodes','off');
close

a=mean(a(1:64,:),2);

shanks=1:8:65;

for i=1:8
    
  REF(i)=find(min(a(shanks(i):(shanks(i+1)-1)))==a );
     REF2(shanks(i):(shanks(i+1)-1))=find(min(a(shanks(i):(shanks(i+1)-1)))==a );
end