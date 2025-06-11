
Signature_of_interest="GABRA1"; % choose signature of interest
Niche_of_interest="CT_region"; % choose region of interest

Threshold_niche=1; %this applies to niche interest
Threshold_Signature=1; %this applies to niche interest




Source_dir="Path to input spreadsheet\"
Source_file="Input_spreadsheet.xlsx"

Target_Dir="Path to output directory\Output"; % Directory for data output
Subfolder="Manders";



Dir=strcat(Target_Dir,"\",Subfolder);




run Import_slide_info_tab.m
clear opts

for tissue_section=1:5; % choose tissues to analyse in input spreadsheet

opts.Sheet = Slide_Info(tissue_section,1)
opts.DataRange = Slide_Info(tissue_section,2)

strcat('tissue',string(tissue_section),'of 5')



run Manders_General_ion.m



%clearvars -except Slide_Info Dir Source_dir Source_file Manders2 Manders_row_names Manders_column_names

%close all


end





