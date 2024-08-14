%laos doer

fn = 'medmod_laos.txt';

strains = ["5","6","7","8","9","10","11","12","13","14","15","16","17","18","19"]; %which strains
ms = [3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3]; %number of harmonics for each

datas = cell(length(strains),1);
for i = 1:length(strains)
    ftr = strcat(fn(1:end-4) , '-strain-' , strains(i));
    datas{i} = RunSPPplus_v2_fnct(ftr,ms(i));    
end
save('laosdata-full.mat','datas');