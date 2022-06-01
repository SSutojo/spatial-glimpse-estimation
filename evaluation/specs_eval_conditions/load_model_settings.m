function [settings_all,param_paths] = load_model_settings(model_settings_txt)
%LOAD_MODEL_SETTINGS reads the text file that contains all parameters and
%arranges them into a cell array is later used to update the parameter struct
%
% All parameters that are not defined in the model_settings_txt are also not
% updated (their values are defined in parameter_settings,
% in the scripts def_prep_params, def_feature_params, etc. )
%
%
% IN
%   model_settings_txt  - name of the txt file as string (e.g. 'conditions_PlotPart2' )
% OUT
%   settings_all        - cell array that contains parameter values that should
%                       be updated for this condition (each column contains a different parameter)
%                       DIMs: number of conditions X number of parameters
%   param_paths         - cell array that defines the struct field of each
%                       updated parameter within the parameter struct.
%                       Corresponds to the columns of settings_all (the
%                       first cell entry param_paths is the struct field of the
%                       parameter in the first column of settings_all)

%% scan txt file, get number of entries, columns, positions of linebreaks etc.
fileID          = fopen([model_settings_txt,'.txt']);
txtcell         = textscan(fileID,'%s');
long_array      = txtcell{1,1};
num_entries     = size(long_array,1);

name_vec_pos        = [];
linebreaks          = [];
feature_column      = [];
rooms_column        = [];
globFeats_column    = [];
for aa = 1:num_entries
    if strcmp(long_array{aa},'%') == 1
        linebreaks      = [linebreaks ; aa];
    end
    if strcmp(long_array{aa},'=') == 1
        name_vec_pos    = aa;
    end
    if strcmp(long_array{aa},'SV_types') == 1  % entries in feature type column are later converted into subcell
        feature_column  = aa-3;
    end
    if strcmp(long_array{aa},'Rooms') == 1  % entries in rooms column are later converted into subcell
        rooms_column  = aa-3;
    end
    if strcmp(long_array{aa},'globFeats') == 1  % entries in rooms column are later converted into subcell
        globFeats_column  = aa-3;
    end
end
num_cols        = linebreaks(2)-4;                      % number of columns
num_rows        = (linebreaks(end-1)-2)/(num_cols+2)-1; % number of rows

%% sort all entries into the new cell array settings_all 
settings_all    = cell(num_rows,num_cols);
for aa = 1:num_rows
    for bb = 1:num_cols
        settings_all{aa,bb} = long_array{4+num_cols + (aa-1)*(num_cols+2) +(bb+1)};
    end
end

%% further change some special types of cell entries
% create subcells in the column with the feature_types/SV_types
for dd = 1:num_rows
    SV_types_long                   = settings_all{dd,feature_column};
    SV_types_split                  = strsplit(SV_types_long,';');
    SV_types                        = cell(size(SV_types_split,2),2);
    for ee = 1:size(SV_types_split,2)
        SV_types(ee,:)              = strsplit(SV_types_split{1,ee},',');
    end
    settings_all{dd,feature_column} = SV_types;
end

if isempty(rooms_column)==0
    for dd = 1:num_rows
        room_types_long                 = settings_all{dd,rooms_column};
        room_types                      = strsplit(room_types_long,',');
        settings_all{dd,rooms_column}   = room_types;
    end
end

if isempty(globFeats_column)==0
    for dd = 1:num_rows
        globFeats_types_long                   = settings_all{dd,globFeats_column};
        globFeats_types_split                  = strsplit(globFeats_types_long,';');
        globFeats_types                        = cell(size(globFeats_types_split,2),2);
        for ee = 1:size(globFeats_types_split,2)
            globFeats_types(ee,:)              = strsplit(globFeats_types_split{1,ee},',');
        end
        settings_all{dd,globFeats_column} = globFeats_types;
    end
end

for ff = 1:num_rows
    for gg = 1:num_cols
        if iscell(settings_all{ff,gg}) == 0  % make sure the entry is no cell
            if ~isnan(str2num(settings_all{ff,gg}))                 % if there is no character in the entry, convert to num
                settings_all{ff,gg} = str2num(settings_all{ff,gg});
            end           
            if strcmp(settings_all{ff,gg},'[]') == 1                % 
                settings_all{ff,gg} = str2num(settings_all{ff,gg});
            end
            if strcmp(settings_all{ff,gg},'true') == 1
                settings_all{ff,gg} = true;
            end
            if strcmp(settings_all{ff,gg},'false') == 1
                settings_all{ff,gg} = false;
            end            
        end
    end
end

%% (field) positions of the changed parameters in the parameter struct
param_long_names    = long_array((name_vec_pos+1):end);  % entries in the second part of the txt file
param_paths         = cell(1,size(param_long_names,1));
for aa = 1:size(param_long_names,1)
    param_paths{aa} = strsplit(param_long_names{aa},'.');
end

