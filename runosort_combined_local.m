function [output] = runosort_combined_local(channel_name)
    
    full_day_path = pwd;
    channels_identified = comb_channels(full_day_path, channel_name);
    
    cd(full_day_path);
    mkdir('session_merged');
    origin = cd('session_merged');
    save('appearances.mat','channels_identified');
    writetable(array2table(channels_identified), 'appearances.csv');
    
    warning('off','MATLAB:MKDIR:DirectoryExists');
    
    for i = 1:size(channels_identified, 1)

                if length(channel_name) == 1
                    out = strcat('channel00', channel_name);
                elseif length(channel_name) == 2
                    out = strcat('channel0', channel_name);
                else
                    out = strcat('channel', channel_name);
                end
                                               
                if channels_identified{i,1} == out
                    split2 = strsplit(channels_identified{i,3}, '/');
                    mkdir(split2{length(split2)});
                    past = cd(split2{length(split2)});
                    mkdir(channels_identified{i,1});
                    cd(channels_identified{i,1});
                    
                    disp('entered');
                    disp(channels_identified{i,1});
                    directories = channels_identified(i,3:2+channels_identified{i,2});
                    
                    data_chunks = cell(1,length(directories));
                    start_indices = cell(2,length(directories));
                    current = 1;
                    
                    for j = 1:length(directories)
                        
                        base_path = cd(directories{j});
                        cd(out);
                        cut = strsplit(directories{j}, '/');

                            start_indices{1,j} = current;
                            start_indices{2,j} = cut{1,length(cut)-1};
                            data_chunks{1,j} = rplhighpass('auto'); % load('rplhighpass.mat');
                            length(data_chunks{1,j}.data.analogData)
                            current = current + length(data_chunks{1,j}.data.analogData);
                            
                        cd(base_path);
                        
                    end
                    
                    disp(current);
                    full_array = zeros(current-1,1);

                    start = 1;
                    for j = 1:length(data_chunks(1,:))

                        disp(size(data_chunks{1,j}.data.analogData));
                        disp(size(full_array(start:start-1+length(data_chunks{1,j}.data.analogData),1)));
                        full_array(start:start-1+length(data_chunks{1,j}.data.analogData),1) = data_chunks{1,j}.data.analogData;
                        start = start + length(data_chunks{1,j}.data.analogData);
                        disp('round');
                    end

                    rw = data_chunks{1,1};

                    rw.data.analogData = single(transpose(full_array));
                %     rw.data = rmfield(rw.data, 'analogTime');
                    rw.data.analogInfo.NumberSamples = current-1;

                    disp(length(rw.data.analogData));

                    save('rplhighpass.mat','rw');
                    save('start_times.mat','start_indices');
                    
                    track_thres = which('runosort.m');
                    
                        fid = fopen(track_thres);

                        tline = fgetl(fid);
                        while ischar(tline)
                            if length(strfind(tline, 'paramsIn.detectionMethod=')) ~= 0
                                if strncmpi(tline, '%', 1) == 0
                                    method_save = tline;
                                end
                            end
                            if length(strfind(tline, 'dp.kernelSize=')) ~= 0
                                if strncmpi(tline, '%', 1) == 0
                                    kernel_save = tline;
                                end
                            end
                            if length(strfind(tline, 'extractionThreshold = ')) ~= 0
                                if strncmpi(tline, '%', 1) == 0
                                    thres_save = tline;
                                end
                            end                            
                            tline = fgetl(fid);
                        end                    
                                       
                    t1 = extractAfter(method_save, '=');
                    t1 = strtrim(extractBefore(t1, ';'))
                    t2 = extractAfter(kernel_save, '=');
                    t2 = strtrim(extractBefore(t2, ';'))
                    t3 = extractAfter(thres_save, '=');
                    t3 = strtrim(extractBefore(t3, ';'))
                    
                    space = ' ';
                    target_clearing = strcat('*detect', t1, 'Thresh', t3, 'kern', t2, '*"');
                    disp(target_clearing);
                    final_command = strcat('for i in `find . -name "', target_clearing, '`; do echo $i; rm -r $i; done;');
                    disp(final_command);
                    unix(char(final_command));
                   
                    
                    RunOSort(pwd);
                    create_pngs(char(strcat('detect', t1, 'Thresh', t3, 'kern', t2)));
                    
                    delete('rplhighpass.mat');
                    
                end
            
    end
end

function [channels_identified] = comb_channels(full_day_path, channel_name)

    origin = cd(full_day_path);
    full_list = dir();
    dirFlags = [full_list.isdir];

    top_folders = full_list(dirFlags);
    top_folders = top_folders(3:length(top_folders));
    
    
    session_names = cell(1,length(top_folders));
    count = 1;
    for i = 1:length(session_names)
        if strncmpi(top_folders(i).name, 'session0', 8) || strncmpi(top_folders(i).name, 'sessioneye', 10)
            session_names{1,count} = top_folders(i).name;
            count = count + 1;
        end
    end
    
    session_names = session_names(1,1:count-1);
        
    day_path = pwd;
    channels_identified = cell(150,count+1);
    identified_count = 0;
    
    for i = 1:length(session_names)
        
        cd(strcat('./',session_names{1,i}));
        
        sub_list = dir();
        dirFlags = [sub_list.isdir];
        sub_folders = sub_list(dirFlags);
        sub_folders = sub_folders(3:length(sub_folders));
        
        array_folders = cell(1,length(sub_folders));
        count = 1;
        
        for j = 1:length(sub_folders)
            if strncmpi(sub_folders(j).name, 'array', 5) == 1
                array_folders{1,count} = sub_folders(j).name;
                count = count + 1;
            end
        end
        array_folders = array_folders(1,1:count-1);
        
        for a = 1:length(array_folders(1,:))
            cd(strcat('./',array_folders{1,a}));
            
            channel_list = dir();
            dirFlags = [channel_list.isdir];
            channel_list = channel_list(dirFlags);
            channel_list = channel_list(3:length(channel_list));

            for k = 1:length(channel_list)
%                 fprintf('%s %s %s\n',session_names{1,i}, array_folders{1,a}, channel_list(k).name);
                if strncmpi(channel_list(k).name, 'channel', 7) == 1
                    found = 0;
                    if identified_count == 0
                        channels_identified{1,1} = channel_list(k).name;
                        channels_identified{1,2} = 1;
                        channels_identified{1,3} = channel_list(k).folder;
                        identified_count = identified_count + 1;
                        found = 1;
                    else
                        for j = 1:identified_count
                            if channels_identified{j,1} == channel_list(k).name
                                channels_identified{j,channels_identified{j,2}+3} = channel_list(k).folder;
                                channels_identified{j,2} = channels_identified{j,2} + 1;
                                found = 1;
                            end
                        end
                    end
                    if found == 0
                        channels_identified{identified_count+1,1} = channel_list(k).name;
                        channels_identified{identified_count+1,2} = 1;
                        channels_identified{identified_count+1,3} = channel_list(k).folder;
                        identified_count = identified_count + 1;
                    end
                end
            end
            
            cd('..');
        end
            
        cd(day_path);
        
    end
    
    channels_identified = channels_identified(1:identified_count,:);
    
%     if string(channel_name) ~= 'all'
%         for i = 1:identified_count
%             if channels_identified{i,1} == string(channel_name)
%                 channels_identified = channels_identified(i,:);
%                 break;
%             end
%         end
%     end
    
    cd(origin);
    
end

