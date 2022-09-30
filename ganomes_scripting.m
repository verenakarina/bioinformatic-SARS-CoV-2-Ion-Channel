path = 'Alpha Genomes'; %change file directory accordingly
T = readtable('dataset/Alpha.csv'); %change file directory accordingly
C = table2cell(T);
num = size(C);
for iteration = 1:num(1)
    titl=char(C{iteration,1});
    my_str = strsplit(titl);
    MN=char(my_str(1));
    if size(MN)>0
        gname=strcat(path, MN);
        fprintf('%d %s \n',iteration,gname)
        getgenbank(MN, 'ToFile', gname)
    end
end