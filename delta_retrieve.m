sout = genbankread('ancestor/SARS-CoV-1-NC_004718.txt');
CDSref=featureparse(sout,'feature','CDS','Sequence',true);


CDStst=6;
sout = genbankread('ancestor/SARS-CoV-2.txt');
CDSq=featureparse(sout,'feature','CDS','Sequence',true);
ts=seqshoworfs(sout.Sequence,'nodisplay', 'true');

OrfStrt=0;
OrfSeq=["ATG"];
 
for m=1:3
    sz=size(ts(m).Stop);
    for k=1:sz(2)
        lns=(ts(m).Stop(k)+2-ts(m).Start(k))+1;
        if (lns >20*3) & (lns<400*3)
            OrfStrt =[OrfStrt, ts(m).Start(k)];
            sq=sout.Sequence( ts(m).Start(k): (ts(m).Stop(k)+2) );
            OrfSeq =[OrfSeq, sq];
        end
    end
end
 
CDSref(CDStst).product
sout=nt2aa(CDSref(CDStst))
queryseq=sout; %(38:133)

[sc,al]=nwalign(queryseq,nt2aa(OrfSeq(2)),'extendgap',1);
scR=-1000.0;
alR=al;
OrfS=0;
sz=size(OrfStrt);
for m=1:(sz(2))
    aa=nt2aa(OrfSeq(m)); 

    [sc,al]=nwalign(queryseq,aa,'extendgap',1);
    if sc>=scR
        alR=al;
        scR=sc;
        OrfS=OrfStrt(m);
    end
end
OrfS
scR
alR

path='Delta_Genomes';
T = readtable('dataset/Delta.csv');
b = table2cell(T);

all_titles=["tmp"];
all_scR=0.0;
all_OrfS=0;
all_ntSequences=["ATG"];
all_indeti=0;
all_simila=0;
all_indel=0;  

numb=size(b)
for ii = 1:numb(1)
    titl=char(b{ii,1});
    my_str = strsplit(titl);
    MN=char(my_str(1));
    if size(MN)>0
        gname=strcat(path, MN);
        sout = genbankread(gname);
        ts=seqshoworfs(sout.Sequence,'nodisplay', 'true'); 
OrfStrt=0;
OrfSeq="ATG"; 
for m=1:3 
    sz=size(ts(m).Stop);
    for k=1:sz(2)
        lns=(ts(m).Stop(k)+2-ts(m).Start(k))+1;
        if (lns >20*3) & (lns<400*3)
            OrfStrt =[OrfStrt, ts(m).Start(k)];
            sq=sout.Sequence( ts(m).Start(k): (ts(m).Stop(k)+2) );
            OrfSeq =[OrfSeq, sq];
        end
    end
end  
[sc,al]=nwalign(queryseq,nt2aa(CDSref(CDStst)),'extendgap',1);
scR=-1000.0;
alR=al;
OrfS=0;
ntSeq=["ATG"];
alignPatterns=[" "]; 
sz=size(OrfStrt);
scores=zeros(1,sz(2));
orflen=zeros(1,sz(2));
for m=1:(sz(2))
    [sc,al]=nwalign(queryseq,  nt2aa(OrfSeq(m)) ,'extendgap',1);
    scores(m) =sc;
    alignPatterns=[alignPatterns, al(2,:)];
    orflen(m)=length(nt2aa(OrfSeq(m)));
    if sc>=scR & OrfStrt(m) >12000
        alR=al;
        scR=sc;
        OrfS=OrfStrt(m);
        ntSeq=OrfSeq(m);
    end
end

seqLen=length(nt2aa(ntSeq));
indeti=length(strfind(alR(2,:),'|'));
simila=length(strfind(alR(2,:),':'));
indel=length(strfind(alR(2,:),' '));
cCount=length(strfind(alR(3,:),'C'));

[scoresS,idx]=sort(scores,'descend');
alPatt = alignPatterns(2:sz(2)+1);

fprintf('%10s, %6.2f, %6.2f, %6.2f, %d, %d, %s, %d, %d, %d,  %d,\n',MN, scoresS(1), scoresS(2), scoresS(3), seqLen, OrfS,alR(2,:),indeti,simila,indel, cCount ); 

alPatt = alignPatterns(2:sz(2)+1);
for m=1:(sz(2))
    if scores(idx(m)) >3.0 & OrfStrt(idx(m)) >22000
    end
end
all_titles=[all_titles, titl];
all_scR=[all_scR,scR];
all_OrfS=[all_OrfS,OrfS]; 
all_ntSequences=[all_ntSequences,ntSeq];
all_indeti=[all_indeti,indeti]; 
all_simila=[all_simila,simila];
all_indel=[all_indel,indel];

    end 
end 

short_all_indeti=0;
short_all_simila=0;
shortlisted=1;

numb=size(all_scR)
for ii = 1:numb(2);
    shortlist=1;
    nb=size(short_all_indeti)
    for k = 1:nb(2)
        if (all_indeti(ii)==short_all_indeti(k) ) && (all_simila(ii)==short_all_simila(k) )
            shortlist=0;
        end
    end
    if shortlist==1
          short_all_indeti=[short_all_indeti,all_indeti(ii)]; 
          short_all_simila=[short_all_simila,all_simila(ii)];        
          shortlisted=[shortlisted,(ii)];        
    end
end

for k = shortlisted 
    ttl=join([all_titles(k), num2str(all_scR(k))]);
    fastawrite('Delta.fas',strcat('Delta',ttl),all_ntSequences(k));
end




