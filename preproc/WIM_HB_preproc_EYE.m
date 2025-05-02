%%
clear
close

%% paths
run(['..' filesep 'localdef_WIM_HB'])
addpath(path_edf2m);

proc_eye = fullfile(wim_preproc, 'eye');
if ~exist(proc_eye, 'dir')
    mkdir(proc_eye)
end

%% Loop on files

files = [dir(['D:\WIM\PBI\sub*' filesep '*.edf']); dir([path_mbi_raweye filesep '*.edf'])];

for n=1:length(files)
    
    fname = files(n).name;
    snum = fname(length('wanderIM_eyelink_s')+(1:3));
    if strcmp(snum, 'MR_')
        sname = 'MWI018';
    else
        sname = ['MWI' snum];
    end

    % import edf into matlab
    fprintf('... converting %s from EDF to .mat file\n', sname)
    myedf = Edf2Mat([files(n).folder filesep fname]);
    
    % clean from useless info
    EL_headers=[];
    EL_headers=myedf.Header;
    EL_headers.Fs=unique(diff(myedf.timeline))*1000'; % in Hertz
    EL_headers.date = [];
    EL_headers.raw = [];

    EL_data=[];
    EL_data.time=myedf.Samples.time;
    EL_data.pupilSize=myedf.Samples.pupilSize;
    EL_data.posX=myedf.Samples.posX;
    EL_data.posY=myedf.Samples.posY;
    
    EL_events=[];
    EL_events.Events.time=myedf.Events.Messages.time;
    EL_events.Events.type=myedf.Events.Messages.info;
    EL_events.StartRec=myedf.Events.Start.time;
    EL_events.EndRec=myedf.Events.End.time;
    
    EL_events.Fix.start=myedf.Events.Efix.start;
    EL_events.Fix.end=myedf.Events.Efix.end;
    EL_events.Fix.duration=myedf.Events.Efix.duration;
    EL_events.Fix.posX=myedf.Events.Efix.posX;
    EL_events.Fix.posY=myedf.Events.Efix.posY;
    EL_events.Fix.pupilSize=myedf.Events.Efix.pupilSize;
    
    EL_events.Blinks.start=myedf.Events.Eblink.start;
    EL_events.Blinks.end=myedf.Events.Eblink.end;
    EL_events.Blinks.duration=myedf.Events.Eblink.duration;
    
    
    EL_events.Sacc.start=myedf.Events.Esacc.start;
    EL_events.Sacc.end=myedf.Events.Esacc.end;
    EL_events.Sacc.duration=myedf.Events.Esacc.duration;
    EL_events.Sacc.posX_start=myedf.Events.Esacc.posX;
    EL_events.Sacc.posY_start=myedf.Events.Esacc.posY;
    EL_events.Sacc.posX_end=myedf.Events.Esacc.posXend;
    EL_events.Sacc.posY_end=myedf.Events.Esacc.posYend;
    EL_events.Sacc.velo=myedf.Events.Esacc.pvel;
    

    % correct doubling
    if length(EL_data.time(1:2:end))==length(EL_data.time(2:2:end)) && min(EL_data.time(1:2:end)-EL_data.time(2:2:end))==0 && max(EL_data.time(1:2:end)-EL_data.time(2:2:end))==0
        warning('time doubled! correcting')
        tim=EL_data.time(1:2:end);
        pupSize=mean([EL_data.pupilSize(1:2:end) EL_data.pupilSize(2:2:end)],2, 'omitnan');
        posX=mean([EL_data.posX(1:2:end) EL_data.posX(2:2:end)],2, 'omitnan');
        posY=mean([EL_data.posY(1:2:end) EL_data.posY(2:2:end)],2, 'omitnan');
        
        EL_data.time=tim;
        EL_data.pupilSize=pupSize;
        EL_data.posX=posX;
        EL_data.posY=posY;
    end
    
    
    % clean ET data
    [data_pupil, filt_data_pupil] = get_EyeLink_cleanpupil(EL_data.pupilSize, EL_headers.Fs, EL_data.time, EL_events);
    EL_data.clean_pupilSize = data_pupil;
    EL_data.filt_pupilSize = filt_data_pupil;

    % save preprocessed data & events
    save([proc_eye filesep sname '_pup'], 'EL_headers','EL_data','EL_events');

end

