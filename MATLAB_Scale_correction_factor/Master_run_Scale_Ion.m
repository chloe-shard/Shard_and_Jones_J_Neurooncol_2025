

%For scale range columns in the non-scaled excel file from A to C
%For other signature analysis using scaled excel file, columns A to n


%Scale_Factor_CytAssist = 4.2839 --> 100um (multiply by 100/4.2839)
%Scale_Factor_Visium = 7.2749 --> 100um (multiply by 100/7.2749)

%signatures to analyze

% # Ref_signature="GABRA1";
% # Neighbourhood_signature="LE_region";


Source_dir="Path to input excel spreadsheet\"
Source_file="Input_spreadsheet.xlsx"

% # Target_Dir="Path to output directory\Data_output"; % Directory for data output
% # Subfolder="GBM_manders_analysis";



% # Dir=strcat(Target_Dir,"\",Subfolder);




run Import_slide_info_tab.m


for tissue_section=1:5; %choose tissues

opts.Sheet = Slide_Info(tissue_section,1)
opts.DataRange = Slide_Info(tissue_section,2)

strcat('tissue',string(tissue_section),'of 1')

run Spot_to_Spot_distance_ion

% # clearvars -except Slide_Info Dir max_k Source_dir Source_file

% #close all

end


figure, hist(r_raw_vector(r_raw_vector<100),1000)
discrete_distance_values=unique(r_raw_vector)'

%figure, hist(discrete_distance_values(discrete_distance_values<100),1000)





