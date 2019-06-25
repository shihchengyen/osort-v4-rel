function [output] = process_channel(current_path)
    
    cd(current_path);
    cut1 = strsplit(current_path, '/');
    channel_name = cut1{1, length(cut1)};
    
    channels_identified = load('../../appearances.mat');
    channels_identified = channels_identified.channels_identified;
    disp(class(channels_identified));
    disp(channels_identified);
    disp(size(channels_identified));
    disp('loaded');
    
    for i = 1:size(channels_identified, 1)
        disp(channels_identified{i,1});
        disp(class(channels_identified{i,1}));
        disp(channel_name);
        disp(class(channel_name));
        if channels_identified{i,1} == channel_name
            directories = channels_identified(i,3:2+channels_identified{i,2});
            break;
        end
    end
    
        data_chunks = cell(1,length(directories));
        start_indices = cell(2,length(directories));
        current = 1;
        
    for j = 1:length(directories)
        command = 'pwd';
        unix(command);
        command = 'scp -P 8398 hippocampus@cortex.nus.edu.sg:/volume1/Hippocampus/Data/picasso-misc/';
        cut = strsplit(directories{j}, '/');
        day_channel_piece = strcat(cut{1,length(cut)-2}, '/', cut{1,length(cut)-1}, '/', cut{1,length(cut)}, '/', channels_identified{i,1});
        command = strcat(command, day_channel_piece, '/rplraw.mat', ' .');
        fprintf('%s\n', command);
        mkdir(cut{1,length(cut)-1});
        cd(cut{1,length(cut)-1});       
        unix(command);
            start_indices{1,j} = current;
            start_indices{2,j} = cut{1,length(cut)-1};
            data_chunks{1,j} = rplraw('auto');
            length(data_chunks{1,j}.data.analogData)
            current = current + length(data_chunks{1,j}.data.analogData);
        cd(current_path);
    end
    
    full_array = zeros(1,current-1);
    
    for j = 1:length(start_indices)-1
        full_array(1,start_indices{1,j}:start_indices{1,j+1}-1) = data_chunks{1,j}.data.analogData;
    end
    
    full_array(1,start_indices{1,length(start_indices)}:current-1) = data_chunks{1,channels_identified{i,2}}.data.analogData;
    rw = data_chunks{1,1};
    
    rw.data.analogData = single(transpose(full_array));
%     rw.data = rmfield(rw.data, 'analogTime');
    rw.data.analogInfo.NumberSamples = current-1;
    
    disp(length(rw.data.analogData));
    
    save('rplraw.mat','rw');
    save('start_times.mat','start_indices');
    
    command = 'rm -r session*';
    unix(command);
    
    RunOSort(pwd);
    create_pngs(pwd);
    
    cd(current_path);
    command = 'rm rplraw.mat';
    unix(command);
    
    command = 'scp -P 8398 -r oSort pngs start_times.mat hippocampus@cortex.nus.edu.sg:/volume1/Hippocampus/Data/picasso-misc/';
    command = strcat(command, cut{1,length(cut)-2}, '/session_merged/', cut{1,length(cut)}, '/', channel_name);
    unix(command);
    
    unix('rm -r oSort');
    unix('rm start_times.mat');
    unix('rm -r pngs');
    
end



