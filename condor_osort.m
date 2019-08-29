function condor_osort()

    disp(dir);
    ls
    directories = dir('rplhighpass_*');
    
        data_chunks = cell(1,length(directories));
        start_indices = cell(2,length(directories));
        current = 1;
        
    for j = 1:length(directories)
        s = ' ';
        unix(['mv', s, directories(j).name, s, 'rplhighpass.mat']);
        sess_id = strsplit(directories(j).name, '_');
        sess_id = sess_id{2};
        sess_id = strsplit(sess_id, '.');
        sess_id = sess_id{1};
            start_indices{1,j} = current;
            start_indices{2,j} = sess_id;
            data_chunks{1,j} = rplhighpass('auto');
            length(data_chunks{1,j}.data.analogData)
            current = current + length(data_chunks{1,j}.data.analogData);
        unix('rm rplhighpass.mat');
    end
    
    disp(current);
    full_array = zeros(current-1,1);
    
    start = 1;
    for i = 1:length(data_chunks(1,:))
        
        disp(size(data_chunks{1,i}.data.analogData));
        disp(size(full_array(start:start-1+length(data_chunks{1,i}.data.analogData),1)));
        full_array(start:start-1+length(data_chunks{1,i}.data.analogData),1) = data_chunks{1,i}.data.analogData;
        start = start + length(data_chunks{1,i}.data.analogData);
        disp('round');
    end
    
    rw = data_chunks{1,1};
    
    rw.data.analogData = single(transpose(full_array));
    rw.data.analogInfo.NumberSamples = current-1;
    
    disp(length(rw.data.analogData));
    
    save('rplhighpass.mat','rw');
    save('start_times.mat','start_indices');
    
    RunOSort(pwd);
    create_pngs(pwd);    
    unix('mkdir pngsfolder');
    unix('mv pngs_* pngsfolder');
    ls
    
    unix('rm rplhighpass.mat');
    
end

