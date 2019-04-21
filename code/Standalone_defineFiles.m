% List the files to sort and their accompanying parameters

dates = {'20181105','20181102','20181101'};

for ii = 1:size(dates,2)
    
    date = dates{ii};

    for jj = 1:123
        if jj<33
            array = 'array01';
        elseif jj>=33 && jj<65
            array = 'array02';
        elseif jj>65 && jj <97
            array = 'array03';
        else 
            array = 'array04';
        end
        channel = num2str(jj);
        channel = horzcat(repmat('0',1,3-length(channel)),channel); 
        
        Standalone_textGUI_demo_new(date,'session01',array,channel);
    end
    
end

% Standalone_textGUI_demo_new('20181105','session01','array01','002');
% Standalone_textGUI_demo_new('20181105','session01','array01','008');
% Standalone_textGUI_demo_new('20181105','session01','array01','024');
% Standalone_textGUI_demo_new('20181105','session01','array01','033');
% Standalone_textGUI_demo_new('20181105','session01','array03','086');
% 
% Standalone_textGUI_demo_new('20181102','session01','array01','010');
% Standalone_textGUI_demo_new('20181102','session01','array01','029');
% Standalone_textGUI_demo_new('20181102','session01','array01','031');
% Standalone_textGUI_demo_new('20181102','session01','array01','032');


