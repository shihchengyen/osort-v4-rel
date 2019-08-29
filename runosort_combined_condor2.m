function [output] = runosort_combined_condor2(channel_name)

    if ~exist('channel_name','var')
        channel_name = 'all';
    end
    
    if isempty(channel_name)
        disp('yes');
        channel_name = 'all';
    end
    
                fid = fopen('~/.bash_profile');

                tline = fgetl(fid);
                while ischar(tline)
                    if strncmpi(tline, 'export GITHUB_MATLAB', 20)
                        tline1 = tline;
                    end
                    if strncmpi(tline, 'export CONDOR_CONFIG', 20)
                        tline2 = tline;
                    end
                    tline = fgetl(fid);
                end

                s = ' ';
                path_1 = strsplit(tline1, '=');
                edit_condor_base = ['cat', s, path_1{2}, '/osort-v4-rel/template_condor_submit.txt', s, '| sed -e ''s/transfer_input_files/transfer_input_files='];
                edit_condor_base = [edit_condor_base, strrep(path_1{2}, '/', '\/'), '\/osort-v4-rel,', s, strrep(path_1{2}, '/', '\/'), '\/DPV,', s, strrep(path_1{2}, '/', '\/'), '\/Hippocampus,' s];
                edit_condor_base = [edit_condor_base, strrep(path_1{2}, '/', '\/'), '\/osort-v4-rel\/setPath_solo.m,', s];
                disp(edit_condor_base);
                path_2 = strsplit(tline2, '=');
                setenv('CONDOR_CONFIG', path_2{2});
    
    
    full_day_path = pwd;
    day = strsplit(full_day_path, '/');
    day = day{1,length(day)};
    channels_identified = comb_channels(full_day_path, channel_name);
    
    cd(full_day_path);
    mkdir('session_merged');
    origin = cd('session_merged');
    save('appearances.mat','channels_identified');
    writetable(array2table(channels_identified), 'appearances.csv');
    
    warning('off','MATLAB:MKDIR:DirectoryExists');
    
    splits = strsplit(string(channel_name), ' ');

    transferred = 0;
    
    for i = 1:size(channels_identified, 1)
        
        if string(channel_name) == 'all'
            
            split2 = strsplit(channels_identified{i,3}, '/');
            mkdir(split2{length(split2)});
            past = cd(split2{length(split2)});
            mkdir(channels_identified{i,1});
            cd(channels_identified{i,1}); % in channel under merged now
            
            generate_txt_command = edit_condor_base;
            first_path = 1;
            for j = 3:2+channels_identified{i,2}
                sess_id = strsplit(channels_identified{i,j}, '/');
                sess_id = sess_id{length(sess_id)-1};
                path_a = strrep(channels_identified{i,j}, '/', '\/');
                
                duplicate_file = cd([channels_identified{i,j}, '/', channels_identified{i,1}]);
                dup_command = 'cp rplhighpass.mat rplhighpass_';
                dup_command = [dup_command, sess_id, '.mat'];
                unix(dup_command);
                cd(duplicate_file);
                
                if first_path == 1
                    generate_txt_command = [generate_txt_command, s, path_a, '\/', channels_identified{i,1}, '\/rplhighpass_', sess_id, '.mat'];
                    first_path = 0;
                else
                    generate_txt_command = [generate_txt_command, ',', s, path_a, '\/', channels_identified{i,1}, '\/rplhighpass_', sess_id, '.mat'];
                end
            end
           
            generate_txt_command = [generate_txt_command, '/g'' > submit_temp.txt'];
            
            disp(generate_txt_command);
            unix(generate_txt_command);
            submit_command = '/usr/local/condor/current/bin/condor_submit submit_temp.txt';
            unix(submit_command);
            unix('rm submit_temp.txt');
            
            cd(past);
            
        else
            
            for j = 1:length(splits(1,:))
                
                if length(splits{1,j}) == 1
                    out = strcat('channel00', splits{1,j});
                elseif length(splits{1,j}) == 2
                    out = strcat('channel0', splits{1,j});
                else
                    out = strcat('channel', splits{1,j});
                end
                
                               
                if channels_identified{i,1} == out
                    split2 = strsplit(channels_identified{i,3}, '/');
                    mkdir(split2{length(split2)});
                    past = cd(split2{length(split2)});
                    mkdir(channels_identified{i,1});
                    cd(channels_identified{i,1});
                    
                    generate_txt_command = edit_condor_base;
                    first_path = 1;
                    for k = 3:2+channels_identified{i,2}
                        sess_id = strsplit(channels_identified{i,k}, '/');
                        sess_id = sess_id{length(sess_id)-1};
                        path_a = strrep(channels_identified{i,k}, '/', '\/');

                        duplicate_file = cd([channels_identified{i,k}, '/', channels_identified{i,1}]);
                        dup_command = 'cp rplhighpass.mat rplhighpass_';
                        dup_command = [dup_command, sess_id, '.mat'];
                        unix(dup_command);
                        cd(duplicate_file);

                        if first_path == 1
                            generate_txt_command = [generate_txt_command, s, path_a, '\/', channels_identified{i,1}, '\/rplhighpass_', sess_id, '.mat'];
                            first_path = 0;
                        else
                            generate_txt_command = [generate_txt_command, ',', s, path_a, '\/', channels_identified{i,1}, '\/rplhighpass_', sess_id, '.mat'];
                        end
                    end
                        
                    generate_txt_command = [generate_txt_command, '/g'' > submit_temp.txt'];

                    disp(generate_txt_command);
                    unix(generate_txt_command);
                    submit_command = '/usr/local/condor/current/bin/condor_submit submit_temp.txt';
                    unix(submit_command);
                    unix('rm submit_temp.txt');

                    cd(past);
                end
            end
        end
    end

    cd(origin);
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

