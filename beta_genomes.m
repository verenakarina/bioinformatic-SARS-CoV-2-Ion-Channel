path = 'Beta Genomes';
T = readtable('dataset/Beta.csv');
b = table2cell(T);
numb=size(b);
for ii = 1:numb(1)
    titl=char(b{ii,1});
    my_str =  strsplit(titl);
    MN=char(my_str(1));
    if size(MN)>0
        gname=strcat(path, MN);
        fprintf('%d %s \n',ii,gname)
        getgenbank(MN, 'ToFile', gname)
    end
end