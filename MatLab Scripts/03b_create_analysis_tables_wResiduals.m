clear
close all

load '/abundance_tables/frequency_table'/KO_freq.mat
load '/abundance_tables/frequency_table'/family_taxa_freq.mat;

load '/Newport'/summary_info/Events_table.mat


%KO
KO_resid = EachAnalysisTable_wResidual(KO_freq,Eventstable);
%Family
Family_resid = EachAnalysisTable_wResidual(family_taxa_freq,Eventstable);

%writetable(KO_resid,'KO_resid.csv')
%writetable(Family_resid,'Family_resid.csv')
