% Copyright (C) 2019 LaureateInstitute for Brain Research
% amayely@laureateinstitute.org
%
% Description: 
% This script runs EEG Preprocessing for EEG Recorded simultaneously with
% fMRI
% For preprocessing steps, please refer to Section 2.6:
% Mayeli, Ahmad, et al. "Real-time EEG artifact correction during fMRI using ICA." Journal of neuroscience methods 274 (2016): 27-37.
% It requires MATLAB and EEGLAB vs. eeglab2019_0
% Authors: Ahmad Mayeli, Kaylee Henry, Chung-Ki Wong
% Contributions: Dr. Jerzy Bodurka's Lab 
% http://www.laureateinstitute.org/jerzy-bodurka.html
% Citation for EEGLAB
%      Delorme A. & Makeig S. (2004), EEGLAB: an open source toolbox for analysis of single-trial EEG dynamics (pdf, 0.7 MB) Journal of Neuroscience Methods 134:9-21
%
% Citation for fmrib_qrsdetect.m and fmrib_pas.m
%   [Niazy06] R.K. Niazy, C.F. Beckmann, G.D. Iannetti, J.M. Brady, and
%   S.M. Smith (2005), Removal of FMRI environment artifacts from EEG data
%   using optimal basis sets. NeuroImage 28 (3), pages 720-737.
%   [Christov04] I.I. Christov (2004), Real time electrocardiogram QRS detection using 
%   combined adaptive threshold, BioMed. Eng. Online 3: 28.
%   [Kim04] K.H. Kim, et al., (2004), Improved ballistocardiac artifact removal 
%   from the electroencephalogram recored in fMRI, J NeouroSience Methods 135: 193-203.

%
% Please Cite:
%   Wong, Chung-Ki, et al. "Automatic cardiac cycle determination directly from EEG-fMRI data by multi-scale peak detection method." Journal of neuroscience methods 304 (2018): 168-184.
%   Mayeli, Ahmad, et al. "Real-time EEG artifact correction during fMRI using ICA." Journal of neuroscience methods 274 (2016): 27-37.
%   Wong, Chung-Ki, et al. "Automatic EEG-assisted retrospective motion correction for fMRI (aE-REMCOR)." Neuroimage 129 (2016): 133.


function APPEAR(varargin)
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
[Filename,indir,infile,outdir,outfile,scntme,dscdt,tr,slpertr,slmkpertr]=deal(varargin{:});
if(ischar(scntme));  scntme = str2double(scntme);  end;
if(ischar(dscdt));  dscdt = str2double(dscdt);  end;
if(ischar(tr));  tr = str2double(tr);  end;
if(ischar(slpertr));  slpertr = str2double(slpertr);  end;
if(ischar(slmkpertr));  slmkpertr = str2double(slmkpertr);  end;
[alphacal,dirMATLAB,outform] = cmdin(Filename);
cdir = pwd;
feeg = 'eeg';
fecg = 'ecg';
fW = 'W';
fEEGpow = 'EEGalpha';
fEEGpowcv = 'EEGalpha_conv';
ficpow = 'ICalpha';
ficpowcv = 'ICalpha_conv';
alphr = [7,13];

if(exist(indir,'dir')~=7)
    error('EEG input directory -- %s -- does not exist\n',indir);
end
if(exist(strcat(indir,'/',infile,'.vhdr'),'file')~=2)
    error('EEG file -- %s/%s.vhdr -- does not exist',indir,infile);
end
if(exist(strcat(indir,'/',infile,'.vmrk'),'file')~=2)
    error('EEG file -- %s/%s.vmrk -- does not exist',indir,infile);
end
if(exist(strcat(indir,'/',infile,'.eeg'),'file')~=2)
    error('EEG file -- %s/%s.eeg -- does not exist',indir,infile);
end
if(exist(outdir,'dir')~=7)
    mkdir(outdir);
end

if(exist(strcat(dirMATLAB,'/MATLAB/ica_linux'),'file')==2)
    system(sprintf('/bin/cp %s/MATLAB/ica_linux %s ; chmod 700 %s/ica_linux',dirMATLAB,cdir,cdir));
else
    error('%s/MATLAB/ica_linux does not exist',dirMATLAB);
end
if(exist(strcat(dirMATLAB,'/MATLAB/eeg_options.m'),'file')==2)
    system(sprintf('/bin/cp %s/MATLAB/eeg_options.m %s ; chmod 700 %s/eeg_options.m',dirMATLAB,cdir,cdir));
else
    error('%s/MATLAB/eeg_options.m does not exist',dirMATLAB);
end
if(exist(strcat(dirMATLAB,'/MATLAB/eeg_optionsbackup.m'),'file')==2)
    system(sprintf('/bin/cp %s/MATLAB/eeg_optionsbackup.m %s ; chmod 700 %s/eeg_optionsbackup.m',dirMATLAB,cdir,cdir));
else
    error('%s/MATLAB/eeg_optionsbackup.m does not exist',dirMATLAB);
end
if(exist(strcat(dirMATLAB,'/MATLAB/cont1.p'),'file')==2)
    system(sprintf('/bin/cp %s/MATLAB/cont1.p %s ; chmod 700 %s/cont1.p',dirMATLAB,cdir,cdir));
else
    error('%s/MATLAB/cont1.p does not exist',dirMATLAB);
end
if(exist(strcat(dirMATLAB,'/MATLAB/conv1.p'),'file')==2)
    system(sprintf('/bin/cp %s/MATLAB/conv1.p %s ; chmod 700 %s/conv1.p',dirMATLAB,cdir,cdir));
else
    error('%s/MATLAB/conv1.p does not exist',dirMATLAB);
end

% load the eeg
[EEG,chlb,mrkA] = EEGLAB(indir,infile,scntme,tr,slmkpertr);
%% Gradient Artifact Correction, DownSampling, and Filtering
% DownSampling, Bandpass filtering, bandstop and sectioning the data
EEG = preproc(EEG,tr,slpertr,mrkA);

% Save the preprocessed data
C = strsplit(infile,'-');
b=sprintf('/%s_FirstPreProc',C{1}');
FirstPreProc_filename=strcat(outdir,b);
%FirstPreProc_filename=sprintf('L:/jbodurka/Ahmad/Evan/MID/Results/%s_FirstPreProc',C{1});
pop_writebva(EEG,FirstPreProc_filename);


%% Try built-in BCG correction on EEGLab (FMRIB plugin)
EEG1=EEG;
EEG2=EEG;
EEG3=EEG;
ECG_chan=32;
[nuECG,deECG] = butter(3,[0.5,15]/(0.5*EEG.srate));
EEG1.data(32,:) = filtfilt(nuECG,deECG,EEG1.data(32,:));
EEG1 = pop_fmrib_qrsdetect(EEG1,ECG_chan,'qrs','no'); % Detect the QRS complexes
evt = EEG1.event;
evtsz = size(evt);
mrkn=0;
%Calculate mean and std of heart rate using fmrib built in function for QRS
%Detection
for ii=1:max(evtsz(1),evtsz(2))
    if(strcmp(evt(ii).type,'qrs'))
        mrkn = mrkn + 1;
        
        QRS(mrkn) = evt(ii).latency;
    end
end
HR1=(250./(diff(QRS)))*60;
stdHR1=std(HR1);
meanHR1=mean(HR1);
QRS_filename=sprintf('%s_QRS.csv',C{1});
csvwrite(QRS_filename,QRS)
ECG=EEG1.data(32,:);
%figure; plot(ECG)
%hold on; plot( QRS, ECG(QRS),  'O');

%% Try CK QRS detection

% determine cardiac period
[pkloc,BCG] = cardperiod(EEG2,slpertr/tr,outfile,outdir,C{1});

EEG2.event=AddMarker(EEG2.event,pkloc,1,0,'R','Response');
EEG2.event=sortlatency(EEG2.event,1);
R_filename=sprintf('%s_R.csv',C{1});
csvwrite(R_filename,pkloc)

%Calculate mean and std of heart rate using CK function for QRS
%Detection
evt = EEG2.event;
evtsz = size(evt);
mrkn=0;
for ii=1:max(evtsz(1),evtsz(2))
    if(strcmp(evt(ii).type,'R'))
        mrkn = mrkn + 1;
        
        Rpeaks(mrkn) = evt(ii).latency;
    end
end
HR2=(250./diff(Rpeaks))*60;
stdHR2=std(HR2);
meanHR2=mean(HR2);


%% Find the more accurate peak detection method
% Find peaks from finger ECG signak
ECG_filename=sprintf('%s_ECG.1D',C{1});
if exist(ECG_filename)
    ECG_Data=csvread(ECG_filename);
    max1=max(ECG_Data(5000:5200));
    max2=max(ECG_Data(8000:8200));
    MinPeakAmp=0.25*min(max1,max2);
    [PKS,PKLOCS]= findpeaks(ECG_Data,'MinPeakDistance',25,'MinPeakHeight',MinPeakAmp);
    PK=round(PKLOCS*250/40);
    HR3=(250./diff(PK))*60;
    stdHR3=std(HR3);
    meanHR3=mean(HR3);
    CK_MeanDif=abs(meanHR3-meanHR2);
    fmrib_MeanDif=abs(meanHR3-meanHR1);
    
    if    CK_MeanDif>(fmrib_MeanDif+0.1)
        EEG1 = pop_fmrib_pas(EEG1,'qrs','mean',21);
        %EEG3 = pop_fmrib_pas(EEG3,'qrs','obs',4);
        EEG=EEG1;
    else
        EEG2 = fmrib_pas(EEG2,pkloc,'mean',21);
        EEG3 = fmrib_pas(EEG3,pkloc,'obs',4);
        EEG=EEG2;
    end
else
    EEG2 = fmrib_pas(EEG2,pkloc,'mean',21);
    EEG3 = fmrib_pas(EEG3,pkloc,'obs',4);
    EEG=EEG2;
end
% eeg-ecg separation
[EEG.data,EEG.nbchan,chlb,ecg1] = eegecgsep(EEG.data,EEG.nbchan,EEG.pnts,chlb);

clear EEG3 EEG1 EEG2
%% Bad Interval Detection
badInter=[];

x = EEG; % save the EEG data to use later for matrix multiplication
SaveTime=EEG.times;
EEG.xmax=EEG.pnts/EEG.srate-1/EEG.srate;
EEG.chanlocs(32)=[];

[EEG, badInter] = pop_rejcont(EEG, 'elecrange',[3:31] ,'freqlimit',[0.5 7] ,'threshold',8,'epochlength',0.5,'contiguous',4,'addlength',0.25,'taper','hamming');


%% Run ICA
n = size(EEG.data,2);
EEG.pnts = n;
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','off');
EEG.chanlocs=loadbvef('BC-MR-32.bvef');
EEG.chanlocs(1)=[]; %removing GRND channel
EEG.chanlocs(1)=[]; %removing REF channel
EEG.chanlocs(32)=[]; %removing ECG channel
figure;
pop_topoplot(EEG, 0, [1:31] ,'',[6 6] ,0,'electrodes','on');

W = EEG.icaweights*EEG.icasphere;
A = inv(W);
% x' = AS' where A is the mixing matrix and x' is the bad interval
% removed data, and S' is calculated below
S_prime = double(W)*double(EEG.data);

% Then we have x, A, S' and we want to solve x=AS for S
% So we do x*W = S
S = double(W)*double(x.data);

% Now we need to get A' which is after we remove the artifactual IC
EEG.times = EEG.times/1000;
EEG.times = [0:0.004:(size(EEG.data,2)-1)*0.004];
%% Classify ICs
[cbicind,saccade,blink,topomap,spectrumy,tpblink,tpsac,smolregion,singchan,muscleloc] = icid(W*double(EEG.data),double(A),double(EEG.data),EEG.srate,EEG.times(end),1);

% Find the columns of A that have artifacts, and that gives us A'
tpblink = logical(tpblink);
tpsac = logical(tpsac);
muscleloc = logical(muscleloc);
singchan = logical(singchan);

% blink or sac
bs = or(tpblink,tpsac);
% muscle or sing chan
ms = or(muscleloc,singchan');
% with bcg
allart = bs+ms+cbicind;
allart = sign(allart);
NeuralICs = find(allart==1);
% A to A'
A(:,NeuralICs) = 0;
A_prime = A;
%% Reconstruct EEG with Inverse ICA
% Multiply A' by S, and that gives us the corrected EEG
EEG.times=SaveTime;
finalEEG=EEG;
finalEEG.times = x.times;
finalEEG.data = A_prime*S; %correct the times
finalEEG.pnts = x.pnts;
finalEEG.event = x.event;
finalEEG.chanlocs = x.chanlocs;
for i = 1:size(badInter,1)
    badlength = badInter(i,2) - badInter(i,1);
    finalEEG.event = AddMarker(finalEEG.event,(badInter(i,1))',badlength,0,'Userdefined','Bad Interval');
end
finalEEG.event=sortlatency(finalEEG.event,1);

% Save final EEG
%corrEEG_filename=sprintf('L:/jbodurka/Ahmad/Evan/MID/Results/%s_corrEEGv10',C{1});
b=sprintf('/%s_corrEEGv10',C{1}');
corrEEG_filename=strcat(outdir,b);
pop_writebva(finalEEG,corrEEG_filename);

end
%% FUNCTIONS:
function output(ind,filename,data,nr,nc)
if(strcmp(ind,'y'))
    fp = fopen(filename,'w');
    for ii = 1:nr
        fprintf(fp,'%10.6f ',data(ii,1:nc));
        fprintf(fp,'\n');
    end
    fclose(fp);
end
end



function outpkt(ind,pkt,filename)
if(strcmp(ind,'y'))
    fp = fopen(filename,'w');
    fprintf(fp,'%10.3f\n',pkt);
    fclose(fp);
end
end



% pwelch calculation
function [psd,frq] = pwelchcal(fn,srate,wsz,ovl)
nr = size(fn,1);
nfft = 2*wsz;
frq = zeros(1,nfft/2+1);
psd = zeros(nr,nfft/2+1);
for ii=1:nr
    [psd(ii,:),frq(:)] = pwelch(fn(ii,:),wsz,ovl,nfft,srate);
end
end



% Independent component analysis
function [A,W] = ICA(Input,rmbad)
tICA = tic;
fprintf('Independent component analysis ...\n');
if(strcmp(rmbad,'y'))
    Input.data = Input.data(1:31,Input.badmot==0);
end
Input = pop_runica(Input, 'icatype', 'runica', 'extended',1,'interrupt','off');
%pop_eegplot(Input,0,1,1)
W = Input.icaweights*Input.icasphere;
%cmps = W * Input.data;
A = inv(W);
toc(tICA);
end



% eeglab input
function [EEG,chlb,mrkA] = EEGLAB(indir,infile,scntme,tr,slmkpertr,slpertr)
% load eeg
fprintf('Loading eeg data from %s/%s\n',indir,infile);
EEG = pop_loadbv(indir,strcat(infile,'.vhdr'));
srate = EEG.srate;
evt = EEG.event;
evtsz = size(evt);
chlb = chlin();

% eeg-fmri data segmentation
smsz = scntme * srate;
tic;
fprintf('Truncating eeg data ...\n');
mrkn = 0;
mrkA = zeros(evtsz(2));
cs = floor(slmkpertr);
frs = slmkpertr - cs;
for ii=evtsz(2):-1:1
    if(strcmp(evt(ii).type,'R128'))
        mrkn = mrkn + 1;
        cs = cs-1;
        mrkA(ii) = evt(ii).latency;
        if(cs==0)
            tr2 = evt(ii).latency - round(tr*srate*frs/slmkpertr);
        end
    end
end
% Added 9/6/19
mrkn = 0;
mrkR128=[];
for ii=1:evtsz(2)
    if(strcmp(evt(ii).type,'R128'))
        mrkn = mrkn + 1;
        
        mrkR128(mrkn) = evt(ii).latency;
    end
end
%  Added
vsTime=1:slmkpertr:mrkn;
for iVol=1:size(vsTime,2)
    mrkAA(iVol)=mrkR128(vsTime(iVol));
end


mrkA = mrkA(mrkA~=0);
tm2 = tr2 + tr*srate - 1;
tm1 = tm2 - smsz + 1;
mrkA = tm1:round(tr*EEG.srate):size(EEG.data,2);
mrkA = mrkA(1:round(scntme/tr));
%{
      if mrkn==slmkpertr*scntme/tr
          mrkA = mrkAA;
      end
%}
EEG.bad = [];
EEG.badmot = [];

toc;
end



% signal preprocessing
function EEG = preproc(EEG,tr,slpertr,mrkA)
dwnsmr = 250;
bpfrq = [1,70];
slfrq = slpertr / tr;
temp = [49, 60, 26, slfrq*(1:1:100)];
temp = temp(temp<=max(2*bpfrq(2),temp(1)));
temp = sort(temp);
bsfrq = zeros(numel(temp),2);
for ii=1:numel(temp)
    bsfrq(ii,1) = temp(ii)-0.5;
    bsfrq(ii,2) = temp(ii)+0.5;
end
nbch = size(EEG.data,1);

% gradient artifact removal
tic;
fprintf('Correcting MRI artifact ...\n');
EEG = fmrib_fastr(EEG,bpfrq(2),1,30,mrkA,0,0,0,numel(mrkA),slpertr,0.03,32,'auto');
EEG.pnts = mrkA(end)+round(tr*EEG.srate)-mrkA(1);
EEG.xmax = EEG.pnts;
EEG.data = EEG.data(:,mrkA(1):mrkA(1)+EEG.pnts-1);
EEG.times = (1:1:EEG.pnts)/EEG.srate;
toc;

% downsampling to 250S/s
tic;
fprintf('Downsampling the eeg data ...\n');
redsmsz = EEG.pnts*dwnsmr/EEG.srate;
sgn1 = zeros(nbch,redsmsz);
for ii=1:nbch
    sgn1(ii,1:redsmsz) = resample(EEG.data(ii,1:EEG.pnts),dwnsmr,EEG.srate,10);
end
EEG.data = sgn1;
RawSampr=EEG.srate;
EEG.srate = dwnsmr;
EEG.pnts = redsmsz;
EEG.xmax = EEG.pnts;
EEG.times = (1:1:EEG.pnts)/EEG.srate;
% mrkA was adjusted for possible triggering time difference, therefore
% it is different from the actual R128 markers.

% Adjusting the marker. Removing R128 markers and adjust the timing after
% segmenting the data and downsampling

%FirstR128=find((mrkA(1)-RawSampr/250)<[EEG.event.latency]& [EEG.event.latency]<(mrkA(1)+RawSampr/250));
First_R128=find((mrkA(1)-800)<[EEG.event.latency]& [EEG.event.latency]<(mrkA(1)+800)); %& strcmp({EEG.event.type},'R128')));

FirstR128=min(First_R128);
if strcmp(EEG.event(FirstR128).type,'R128')
    FirstR128=FirstR128;
else
    FirstR128=First_R128(2);
end
str=EEG.event(FirstR128(1)).latency;
EEG.event=EEG.event(FirstR128(1):size(EEG.event,2));

CurrEventSize=size(EEG.event,2);

EEG.event=AddMarker(EEG.event,mrkA,1,0,'VS','Response');
for i=1:size(EEG.event,2)
    EEG.event(i).latency= round((EEG.event(i).latency-str)/(RawSampr/EEG.srate));
    EEG.event(i).bvmknum= round((EEG.event(i).bvmknum-FirstR128+1));
    EEG.event(i).urevent= round((EEG.event(i).urevent-FirstR128+1));
end
cnt=0;
for i=1:size(EEG.event,2)
    if strcmp(EEG.event(i).type,'R128')==1 || strcmp(EEG.event(i).code,'SyncStatus')==1
        cnt=cnt+1;
        DltMark(cnt)=i;
    end
    
end
EEG.event(DltMark)=[];

EEG.event=sortlatency(EEG.event,1);
EEG.urevent=EEG.event;
EEG.urevent=rmfield(EEG.urevent,'urevent');
clear ReducedSignal SrtEventTable NumEvent EventTable cnt i CurrEventSize;
toc;
%%
% bandpass filtering
tic;
fprintf('Band-pass frequency filtering the eeg data ...\n');
bpfrq(2) = min(bpfrq(2),0.5*EEG.srate);
EEG.data = eegfilt(EEG.data,EEG.srate,bpfrq(1),0,0,3*fix(EEG.srate/1),0,'fir1',0);
EEG.data = eegfilt(EEG.data,EEG.srate,0,bpfrq(2),0,3*fix(EEG.srate/1),0,'fir1',0);
toc;

% bandstop filtering at multiple acquisition frequencies, 26Hz, 60Hz
tic;
bsfrq = bsfrq(bsfrq(:,2)<0.5*EEG.srate,1:2);
bsfrq2 = size(bsfrq);
fprintf('Band-stop frequency filtering the eeg data ...\n');
for iibs = 1:bsfrq2(1)
    Wn = bsfrq(iibs,1:2)/(0.5*EEG.srate);
    [nu,de] = butter(3,Wn,'stop');
    for ii=1:nbch
        EEG.data(ii,:) = filtfilt(nu,de,EEG.data(ii,:));
    end
end
toc;
end



function mbpfn = icsel(fcn,srate,cbind,acqfrq)
[nu,de] = butter(3,[2,2*acqfrq]/srate);
for ii=find(cbind==1)
    fcn(ii,:) = filtfilt(nu,de,fcn(ii,:));
end
sgn = zeros(1,size(fcn,1));
cnt = zeros(1,size(fcn,1));
lccnt = zeros(size(fcn,1),floor(size(fcn,2)/floor(4*srate)));
mn1 = mean(fcn,2);
std1 = std(fcn,[],2);
pssd = zeros(1,size(fcn,1));
ngsd = zeros(1,size(fcn,1));
for ii = find(cbind==1)
    fcn(ii,:) = fcn(ii,:)/std1(ii);
    mn1(ii) = mean(fcn(ii,:));
    std1(ii) = std(fcn(ii,:));
    pssd(ii) = mean( fcn(ii, fcn(ii,:)>=mn1(ii)+std1(ii) & fcn(ii,:)<=mn1(ii)+4*std1(ii) ) );
    ngsd(ii) = mean( fcn(ii, fcn(ii,:)<=mn1(ii)-std1(ii) & fcn(ii,:)>=mn1(ii)-4*std1(ii) ) );
    if(abs(pssd(ii))>=abs(ngsd(ii)))
        sgn(ii) = 1;
    elseif(abs(pssd(ii))<abs(ngsd(ii)))
        sgn(ii) = -1;
    end
    fcn(ii,:) = sgn(ii)*fcn(ii,:);
    mn1(ii) = sgn(ii)*mn1(ii);
    Fn = fcn(ii,:);  Fn(Fn>mn1(ii)+4*std1(ii))=mn1(ii)+4*std1(ii);
    for jj=1:floor(size(fcn,2)/floor(4*srate))
        seg = Fn((jj-1)*floor(4*srate)+1:jj*floor(4*srate));
        segpk = findpeaks(seg);
        segpk = sort(segpk,'descend');
        segval = sort(seg,'descend');
        lccnt(ii,jj) = (mean(segpk(1:round(0.1*numel(segpk))))-mean(segval(round(0.1*numel(segval)):end)))/std(segval(round(0.1*numel(segval)):end));
    end
    cnt(ii) = mean(lccnt(ii,:));
end
ind1 = find(cbind==1);
[~,ind] = max(cnt(ind1));
mbpfn = fcn(ind1(ind),:);
end



% cardiac cycle determination
function [pkloc, mbpfn] = cardperiod(input1,acqfrq,sess,outdir,EXP)
tic;
fprintf('Determining cardiac period ...\n');
input1.data = input1.data(1:31,:);
[A,W] = ICA(input1,'n');
cbind = icid(W*input1.data,A,input1.data,input1.srate,input1.times(end),0);
mbpfn = icsel(W*input1.data,input1.srate,cbind,acqfrq);
BCG_filename=sprintf('%s_BCG.csv',EXP);
csvwrite(BCG_filename,mbpfn)
pkloc = mbp(mbpfn,input1.srate,acqfrq);
toc;
end



% data segmentation
function [EEG,pkloc] = segmnt(EEG,pkloc,dscdt)
EEG.data = EEG.data(:,dscdt*EEG.srate+1:EEG.xmax);
EEG.bad = EEG.bad(:,dscdt*EEG.srate+1:EEG.xmax);
EEG.badmot = EEG.badmot(dscdt*EEG.srate+1:EEG.xmax);
pkloc = pkloc - round(dscdt*EEG.srate);
pkloc = pkloc(pkloc>0);
EEG.pnts = size(EEG.data,2);
EEG.xmax = size(EEG.data,2);
EEG.times = (1:1:EEG.xmax)/EEG.srate;
cnt=0;
for i=1:max(size(EEG.event,1),size(EEG.event,2))
    EEG.event(i).latency=EEG.event(i).latency-dscdt*EEG.srate;
    if EEG.event(i).latency<0
        cnt=cnt+1;
        DltMark(cnt)=i;
    end
end
if cnt>0
    EEG.event(DltMark)=[];
end
end



% eeg-ecg signal separation
function [sgn,nbch,chlb,sgn1] = eegecgsep(sgn,nbch,smsz,chlb)
sgn1 = sgn(nbch,1:smsz);
sgn = sgn(1:nbch-1,1:smsz);
tmp = cell(nbch-1,3);
for ii=1:nbch-1
    tmp{ii,1} = chlb{ii,1};
end
chlb = tmp;
nbch = nbch-1;
end



% bad segment definition
function badsegind = badseg(EEG,chlb)
[EEG.data,EEG.nbchan,~,~] = eegecgsep(EEG.data,EEG.nbchan,EEG.pnts,chlb);
ext = 0.3;
dscdt = 0.3;
sgn = EEG.data;
sgn = eegfilt(sgn,EEG.srate,1,0,0,3*fix(EEG.srate),0,'fir1',0);
nbch = EEG.nbchan;
smsz = EEG.pnts;
sgnmn = mean(EEG.data,2);
sgnstd = std(EEG.data,0,2);
durpt = floor(0.04*EEG.srate);
extpt = floor(ext*EEG.srate);
dscdpt = floor(dscdt*EEG.srate);
badsegind = zeros(size(sgn));
for ii=1:nbch
    ind = ( sgn(ii,:)>sgnmn(ii)+4*sgnstd(ii) | sgn(ii,:)<sgnmn(ii)-4*sgnstd(ii) );
    [sctn,sctr] = cont1(ind);
    if(sctn>0)
        for jj=1:sctn
            if(sctr(jj,2)-sctr(jj,1)>=durpt && max(abs(sgn(ii,sctr(jj,1):sctr(jj,2))))>200)
                sctr(jj,1) = max(1,sctr(jj,1) - extpt);
                sctr(jj,2) = min(smsz,sctr(jj,2) + extpt);
                ind(sctr(jj,1):sctr(jj,2)) = 1;
            else
                ind(sctr(jj,1):sctr(jj,2)) = 0;
            end
        end
        [sctn,sctr] = cont1(ind);
        if(sctn>0)
            for jj=1:sctn
                ltpt = sctr(jj,1);
                while( ltpt>1 && sgn(ii,ltpt)*sgn(ii,ltpt-1)>0 )
                    ind(ltpt-1) = 1;
                    ltpt = ltpt - 1;
                end
                rtpt = sctr(jj,2);
                while( rtpt<smsz && sgn(ii,rtpt)*sgn(ii,rtpt+1)>0 )
                    ind(rtpt+1) = 1;
                    rtpt = rtpt + 1;
                end
            end
            badsegind(ii,:) = ind;
        end
    end
end

for ii=1:nbch
    [sctn,sctr] = cont1(~badsegind(ii,:));
    for jj=1:sctn
        if(numel(sctr(jj,1):sctr(jj,2))<dscdpt)
            badsegind(ii,sctr(jj,1):sctr(jj,2)) = 1;
        end
    end
end
end



% suspected motion segment
function bd2ind = badsegmo(EEG)
sdeind = ones(1,EEG.nbchan);
sdeind([17:19,21:24,32]) = 0;
Index = EEG.bad(find(sdeind==1),:);
Index(Index(:, sum(Index(:,:),1) < 2)==1) = 0;
bd2ind = ( sum(Index,1)>=6 );
end



% downsampling
function dwndata = dwnsamp(data,orig,new)
dwndata = zeros(size(data,1)*orig/new,size(data,2));
for ii=1:size(data,2)
    dwndata(:,ii) = resample(data(:,ii),1,new/orig);
end
end



% electrode labels
function [chlb] = chlin()
chlb = {'Fp1';'Fp2';'F3';'F4';'C3';'C4';'P3';'P4';...
    'O1';'O2';'F7';'F8';'T7';'T8';'P7';'P8';...
    'Fz';'Cz';'Pz';'Oz';'FC1';'FC2';'CP1';'CP2';...
    'FC5';'FC6';'CP5';'CP6';'TP9';'TP10';'POz';'ECG'};
end



% import parameters
function [alphacal,dirMATLAB,outform] = cmdin(filename)
if(exist(filename,'file')~=2)
    error(sprintf('%s does not exist',filename));
end
fp = fopen(filename);
formatSpec = 'AlphaPower = %s MATLABdir = %s OutputFormat = %s';
C = textscan(fp,formatSpec,'Delimiter','\n','CollectOutput',true);
fclose(fp);
alphacal = C{1}{1};
dirMATLAB = C{1}{2};
outform = C{1}{3};
if(exist(strcat(dirMATLAB,'/MATLAB'),'dir')~=7)
    error(sprintf('%s/MATLAB does not exist',dirMATLAB));
end
end



function [dst] = lne2pt(lnpt1, lnpt2, pt)
A = lnpt2 - lnpt1;
r = lnpt1 - pt;
rxA = cross(r,A);
magA = sqrt(sum(A.*A));
dst = sqrt(sum(rxA.*rxA))/magA;
end



% bounded linear interpolation
function [yy0] = linintep(xx,yy,xx0)
yy0 = (yy(1)-yy(2))/(xx(1)-xx(2))*xx0 + (yy(2)*xx(1)-yy(1)*xx(2))/(xx(1)-xx(2));
if(yy0>max(yy))
    yy0 = max(yy);
elseif(yy0<min(yy))
    yy0 = min(yy);
end
end


% calculate few statistics
function [fave,fstd,fvar,fskw,fkurt] = stat(fn)
szfn = size(fn);
nr = szfn(1);
nc = szfn(2);
fave = zeros(1,nr);
fstd = zeros(1,nr);
fvar = zeros(1,nr);
fskw = zeros(1,nr);
fkurt = zeros(1,nr);
fave(:) = mean(fn(:,1:nc),2);
fstd(:) = std(fn(:,1:nc),0,2);
fvar(:) = var(fn(:,1:nc),0,2);
for ii=1:nr
    fskw(ii) = sum((fn(ii,1:nc)-fave(ii)).^3)/(nc*fvar(ii)*fstd(ii));
    fkurt(ii) = sum((fn(ii,1:nc)-fave(ii)).^4)/(nc*fvar(ii)*fvar(ii))-3;
end
end



% function max/min
function [pkl,lftmin,rgtmin,pkn,dpt,dptave,dptstd] = fninfo(fn)
[minv,minlc] = findpeaks(-fn);
nbmin = size(minv);
pkn = nbmin(2)-1;

pkl = zeros(1,pkn);
lftmin = zeros(1,pkn);
rgtmin = zeros(1,pkn);
lftmin(1:pkn) = minlc(1:pkn);
rgtmin(1:pkn) = minlc(2:pkn+1);
for ii = 1:pkn
    [~,pkl(ii)] = max(fn(lftmin(ii):rgtmin(ii)));
    pkl(ii) = pkl(ii) + lftmin(ii) - 1;
end
dpt = fn(pkl)-0.5*(fn(lftmin)+fn(rgtmin));

dptave = mean(dpt);
dptstd = std(dpt);
end



% automatic cardiac cycle determination directly from EEG-fMRI data
function pkloc = mbp(Fn,srate,acqfrq)
tmbp = tic;

% estimation of cardiac cycle
spt = abs(fft(Fn.*(Fn>=0))).^2;
sgnpt = fft(spt);
tmwin = zeros(size(sgnpt));
hannwin = hann(2*round(8*srate),'periodic');
tmwin(1:round(8*srate)) = hannwin(round(8*srate)+1:2*round(8*srate));
tmwin(length(tmwin)-round(8*srate)+1:length(tmwin)) = hannwin(1:round(8*srate));
sgnpt = sgnpt.*tmwin;
smspt = real(ifft(sgnpt));
smspt = smspt/max(smspt(round(0.5*size(Fn,2)/srate):round(8*size(Fn,2)/srate)));
spt = spt/max(spt(round(0.5*size(Fn,2)/srate):round(8*size(Fn,2)/srate)));
[~, allpklocs] = findpeaks(smspt);
[~, allvalocs] = findpeaks(-smspt);
allvalocs(allvalocs*srate/size(Fn,2)>8) = [];
allpklocs(allpklocs*srate/size(Fn,2)>8) = [];
if(allpklocs(1)<allvalocs(1))
    allpklocs(1) = [];
end
if(allpklocs(end)>allvalocs(end))
    allpklocs(end) = [];
end
allpklocs = zeros(1,numel(allvalocs)-1);
for ii=1:numel(allvalocs)-1
    [~,ind] = max(smspt(allvalocs(ii):allvalocs(ii+1)));
    allpklocs(ii) = ind + allvalocs(ii) - 1;
end

%  find fundamental, bp1, first harmonic and reference frequencies
pkloc = [];
hbr = 120;
while(numel(pkloc)==0 && hbr<=150)
    pkloc = allpklocs(allpklocs>=floor(0.5/srate*length(Fn))&allpklocs<=floor((hbr/60)/srate*length(Fn)));
    hbr = hbr + 7.5;
end
selval = 0;
for jj=1:numel(pkloc)
    valoc = allvalocs(allvalocs<pkloc(jj));
    [~, index] = min(abs(valoc - pkloc(jj)));
    lfvaloc = valoc(index);
    valoc = allvalocs(allvalocs>pkloc(jj));
    [~, index] = min(abs(valoc - pkloc(jj)));
    rgvaloc = valoc(index);
    data = sort(spt(lfvaloc:rgvaloc),'descend');
    tmpv = smspt(pkloc(jj));
    if(tmpv>selval)
        selval = tmpv;
        sptpkloc = pkloc(jj);
        smsplfvaloc = lfvaloc;
        smsprgvaloc = rgvaloc;
    end
end

valoc = allvalocs(allvalocs<=smsplfvaloc-1);
pkloc = allpklocs(allpklocs<=sptpkloc-1);
if(numel(valoc)>0 && numel(pkloc)>0)
    valoc = sort(valoc,'descend');
    pkloc = sort(pkloc,'descend');
    val1 = smspt(sptpkloc);
    val2 = smspt(smsplfvaloc);
    jj = 1;
    while(and(jj<=numel(valoc),jj<=numel(pkloc)) && smspt(pkloc(jj))>=0.15*smspt(sptpkloc) && smspt(pkloc(jj))<val1 && smspt(valoc(jj))<val2 && pkloc(jj)>round(0.5*size(Fn,2)/srate))
        val1 = smspt(pkloc(jj));
        val2 = smspt(valoc(jj));
        smsplfvaloc = valoc(jj);
        jj = jj + 1;
    end
end
rg1 = max(sptpkloc*2-0.5*(smsprgvaloc-smsplfvaloc),smsprgvaloc);
rg2 = min(sptpkloc*2+0.5*(smsprgvaloc-smsplfvaloc),2.5*sptpkloc);
pkloc = allpklocs( allpklocs >= rg1 & allpklocs <= rg2 );
while(numel(pkloc)==0)
    rg1 = max(rg1 - 0.05*(smsprgvaloc-smsplfvaloc),smsprgvaloc);
    rg2 = min(rg2 + 0.05*(smsprgvaloc-smsplfvaloc),2.5*sptpkloc);
    pkloc = allpklocs( allpklocs >= rg1 & allpklocs <= rg2 );
end
[~,ind] = max(smspt(pkloc));
hmpkloc = pkloc(ind);
hmlfvaloc = max(allvalocs(allvalocs<hmpkloc));
valoc = allvalocs(allvalocs>=smsprgvaloc+1 & allvalocs<hmpkloc);
pkloc = allpklocs(allpklocs>=sptpkloc+1 & allpklocs<hmpkloc);
val1 = smspt(sptpkloc);
val2 = smspt(smsprgvaloc);
jj = 1;
while( and(jj<=numel(valoc),jj<=numel(pkloc)) && smspt(pkloc(jj))>=0.15*smspt(sptpkloc) && smspt(pkloc(jj))<val1 && smspt(valoc(jj))<val2 && valoc(jj)<=hmlfvaloc )
    val1 = smspt(pkloc(jj));
    val2 = smspt(valoc(jj));
    smsprgvaloc = valoc(jj);
    rg1 = max(sptpkloc*2-0.5*(smsprgvaloc-smsplfvaloc),smsprgvaloc);
    rg2 = min(sptpkloc*2+0.5*(smsprgvaloc-smsplfvaloc),2.5*sptpkloc);
    pkloc = allpklocs( allpklocs >= rg1 & allpklocs <= rg2 );
    while(numel(pkloc)==0)
        rg1 = max(rg1 - 0.05*(smsprgvaloc-smsplfvaloc),smsprgvaloc);
        rg2 = min(rg2 + 0.05*(smsprgvaloc-smsplfvaloc),2.5*sptpkloc);
        pkloc = allpklocs( allpklocs >= rg1 & allpklocs <= rg2 );
    end
    [~,ind] = max(smspt(pkloc));
    hmpkloc = pkloc(ind);
    hmlfvaloc = max(allvalocs(allvalocs<hmpkloc));
    jj = jj + 1;
end
bp1bpfrq = [smsplfvaloc smsprgvaloc] * srate/size(Fn,2);

lmt = 0.5*size(Fn,2)/srate;
loc1 = sptpkloc;
val = mean(spt(sptpkloc-round(0.05*size(Fn,2)/srate):sptpkloc+round(0.05*size(Fn,2)/srate)));
jj = -1;
rge = [sptpkloc+round((jj*0.05-0.05)*size(Fn,2)/srate),sptpkloc+round((jj*0.05+0.05)*size(Fn,2)/srate)];
while(rge(1)>lmt)
    if(val>mean(spt(rge(1):rge(2)))&&mean(rge)<=smsplfvaloc)
        loc1 = mean(rge);
        val = mean(spt(rge(1):rge(2)));
    end
    jj = jj - 1;
    rge = [sptpkloc+round((jj*0.05-0.05)*size(Fn,2)/srate),sptpkloc+round((jj*0.05+0.05)*size(Fn,2)/srate)];
end
loc1 = min(loc1,smsplfvaloc);
lmt = hmpkloc-round((smsprgvaloc-smsplfvaloc)/2);
loc2 = sptpkloc;
val = mean(spt(sptpkloc-round(0.05*size(Fn,2)/srate):sptpkloc+round(0.05*size(Fn,2)/srate)));
jj = 1;
rge = [sptpkloc+round((jj*0.05-0.05)*size(Fn,2)/srate),sptpkloc+round((jj*0.05+0.05)*size(Fn,2)/srate)];
while(rge(2)<lmt)
    if(val>mean(spt(rge(1):rge(2)))&&mean(rge)>=smsprgvaloc)
        loc2 = mean(rge);
        val = mean(spt(rge(1):rge(2)));
    end
    jj = jj + 1;
    rge = [smsprgvaloc+round((jj*0.05-0.05)*size(Fn,2)/srate),smsprgvaloc+round((jj*0.05+0.05)*size(Fn,2)/srate)];
end
loc2 = max(loc2,smsprgvaloc);
refloc3 = loc2;
refloc1 = loc1;
refloc2 = loc1+0.5*(sptpkloc-loc1);

% calculate bp1
Fnp = Fn;  Fnp(Fnp<0) = 0;  Fnp = Fnp.*abs(Fnp);
ftbp1 = fft(Fnp);
ftbp1(1:floor(bp1bpfrq(1)*size(Fn,2)/srate)) = 0;
ftbp1(length(Fn)-floor(bp1bpfrq(1)*size(Fn,2)/srate)+1:length(Fn)) = 0;
ftbp1(floor(bp1bpfrq(2)*size(Fn,2)/srate)+1:length(Fn)-floor(bp1bpfrq(2)*size(Fn,2)/srate)) = 0;
bp1 = real(ifft(ftbp1));
bp1 = bp1/std(bp1);
[bp1pkloc,bp1lfloc,bp1rgloc,~,~,~,~] = fninfo(bp1);
bp1lfbd = bp1lfloc;
bp1rgbd = bp1rgloc;
nbpk = size(bp1pkloc,2);
ii = 2;
while(1)
    if(ii>=nbpk)
        break;
    end
    prvpk = bp1pkloc(max(ii-10,1):ii-1);
    prvlfloc = bp1lfloc(max(ii-10,1):ii-1);
    prvrgloc = bp1rgloc(max(ii-10,1):ii-1);
    dptmn = mean(bp1(prvpk)-0.5*(bp1(prvlfloc)+bp1(prvrgloc)));
    cond11a = ( bp1(bp1pkloc(ii))-bp1(bp1lfloc(ii))<0.2*dptmn && bp1pkloc(ii)-bp1pkloc(ii-1)<1/(refloc3/size(Fn,2)) );
    cond11b = ( bp1(bp1pkloc(ii))-bp1(bp1rgloc(ii))<0.2*dptmn && bp1pkloc(ii+1)-bp1pkloc(ii)<1/(refloc3/size(Fn,2)) );
    cond11c = ( or(bp1(bp1pkloc(ii))-bp1(bp1lfloc(ii))<0.2*dptmn,bp1(bp1pkloc(ii))-bp1(bp1rgloc(ii))<0.2*dptmn) && bp1rgloc(ii)-bp1lfloc(ii)<1/(refloc3/size(Fn,2)) );
    cond12 = ( bp1pkloc(ii+1)-bp1pkloc(ii-1)<srate/bp1bpfrq(1) );
    if( (cond11a && cond12) || (cond11b && cond12) || (cond11c && cond12) )
        if( cond11a )
            bp1rgbd(ii-1) = bp1pkloc(ii);
            bp1lfbd(ii+1) = bp1pkloc(ii);
            bp1rgloc(ii-1) = bp1rgloc(ii);
        elseif( cond11b )
            bp1lfbd(ii+1) = bp1pkloc(ii);
            bp1rgbd(ii-1) = bp1pkloc(ii);
            bp1lfloc(ii+1) = bp1lfloc(ii);
        elseif( cond11c)
            bp1rgbd(ii-1) = bp1pkloc(ii);
            bp1lfbd(ii+1) = bp1pkloc(ii);
        end
        bp1pkloc(ii) = [];
        bp1lfloc(ii) = [];
        bp1rgloc(ii) = [];
        bp1lfbd(ii) = [];
        bp1rgbd(ii) = [];
        nbpk = nbpk - 1;
        ii = ii - 1;
    end
    ii = ii + 1;
end

%  calculate bp2
[nu,de] = butter(3,[1,acqfrq]/(0.5*srate));
bp2 = filtfilt(nu,de,Fn);
for ii = 2:numel(bp1pkloc)-1
    rge = [bp1lfbd(ii),bp1rgbd(ii)-1];
    bp2(rge(1):rge(2)) = bp2(rge(1):rge(2))/std(bp2(rge(1):rge(2)));
end
rge = [1,bp1lfbd(2)-1];
bp2(rge(1):rge(2)) = bp2(rge(1):rge(2))/std(bp2(rge(1):rge(2)));
rge = [bp1rgbd(numel(bp1pkloc)-1),size(bp2,2)];
bp2(rge(1):rge(2)) = bp2(rge(1):rge(2))/std(bp2(rge(1):rge(2)));
[bp2pkloc,~,~,~,~,~,~] = fninfo(bp2);
bp2pkamp = bp2(bp2pkloc);
bp2lfslp = zeros(size(bp2pkloc));
bp2rgslp = zeros(size(bp2pkloc));
for ii=1:numel(bp2pkloc)
    bp2lfslp(ii) = bp2pkamp(ii)-min(bp2(max(1,bp2pkloc(ii)-round(0.05*srate)):bp2pkloc(ii)));
    bp2rgslp(ii) = bp2pkamp(ii)-min(bp2(bp2pkloc(ii):min(bp2pkloc(ii)+round(0.05*srate),size(bp2,2))));
end
Bp2 = struct('data',bp2,'pkloc',bp2pkloc,'pkamp',bp2pkamp,'lfslp',bp2lfslp,'rgslp',bp2rgslp,'selpkloc',zeros(1,numel(bp1pkloc)));

% select one peak in each estimated cycle
for ii=1:numel(bp1pkloc)
    itvl = [bp1lfbd(ii),bp1rgbd(ii)];
    rge = (max(ii-10,1):ii-1);
    coritvl = zeros(numel(rge),2);
    coritvlpkloc = zeros(1,numel(rge));
    for jj=1:numel(rge)
        coritvl(jj,:) = [bp1lfbd(rge(jj)),bp1rgbd(rge(jj))];
        coritvlpkloc(jj) = Bp2.selpkloc(Bp2.selpkloc>=coritvl(jj,1) & Bp2.selpkloc<coritvl(jj,2));
    end
    [pkind,sel,selflg,Bp2] = selpk2(Bp2,itvl,coritvl,coritvlpkloc,[mean(Bp2.data)+2*std(Bp2.data),mean(max([Bp2.lfslp;Bp2.rgslp],[],1))+2*std(max([Bp2.lfslp;Bp2.rgslp],[],1))],[1 1 1 1]);
    Bp2.selpkloc(ii) = Bp2.pkloc(pkind(find(sel==1)));
end

% add 1st and last peak if needed
nghb = [numel(bp1pkloc),numel(bp1pkloc)-10; 1,10];
bd = [bp1rgbd(end),size(Fn,2);bp1lfbd(1),1];
for ii = 1:size(nghb,1)
    if(max(bd(ii,:))-min(bd(ii,:))>3)
        [~, bp1pkloc1] = findpeaks(bp1(min(bd(ii,:)):max(bd(ii,:))));
        if(numel(bp1pkloc1)>0)
            bp1pkloc1 = bp1pkloc1 + min(bd(ii,:)) - 1;
            [~,ind] = max(bp1(bp1pkloc1));
            bp1pkloc1 = bp1pkloc1(ind);
            Bp1lfrs = bp1(bp1pkloc)-bp1(bp1lfloc);
            Bp1rgrs = bp1(bp1pkloc)-bp1(bp1rgloc);
            dpt = 0.5*(Bp1lfrs(min(nghb(ii,:)):max(nghb(ii,:)))+Bp1rgrs(min(nghb(ii,:)):max(nghb(ii,:))));
            if(1)
                itvl = [min(bd(ii,:)),max(bd(ii,:))];
                rge = nghb(ii,1:2);
                coritvl = zeros(numel(rge),2);
                coritvlpkloc = zeros(1,numel(rge));
                for jj=1:numel(rge)
                    coritvl(jj,:) = [bp1lfbd(rge(jj)),bp1rgbd(rge(jj))];
                    coritvlpkloc(jj) = Bp2.selpkloc(Bp2.selpkloc>=coritvl(jj,1) & Bp2.selpkloc<coritvl(jj,2));
                end
                [pkind,sel,selflg,Bp2] = selpk2(Bp2,itvl,coritvl,coritvlpkloc,[mean(Bp2.data)+2*std(Bp2.data),mean(max([Bp2.lfslp;Bp2.rgslp],[],1))+2*std(max([Bp2.lfslp;Bp2.rgslp],[],1))],[1 1 1 1]);
                if(bd(ii,2)==size(Fn,2) && all(selflg(find(sel==1),:)==1))
                    Bp2.selpkloc = [Bp2.selpkloc,Bp2.pkloc(pkind(find(sel==1)))];
                elseif(bd(ii,2)==1 && all(selflg(find(sel==1),:)==1))
                    Bp2.selpkloc = [Bp2.pkloc(pkind(find(sel==1))),Bp2.selpkloc];
                end
            end
        end
    end
end

% cycle adjustment
sep0 = 0.5*size(Fn,2)/sptpkloc;
sep1 = size(Fn,2)/refloc3;
sep2 = size(Fn,2)/refloc2;
sep3 = size(Fn,2)/refloc1;
ii=2;
flag = 1;
while(flag)
    ii = ii + 1;
    if(ii>=3 && ii<=numel(Bp2.selpkloc)-1)
        dis1 = Bp2.selpkloc(ii)-Bp2.selpkloc(ii-1);
        dis0 = Bp2.selpkloc(ii-1)-Bp2.selpkloc(ii-2);
        dis2 = Bp2.selpkloc(ii+1)-Bp2.selpkloc(ii);
        cond21a = ( dis1<sep1 && or(dis0>sep2,dis2>sep2) && (dis0+dis1+dis2)>=2.5*size(Fn,2)/sptpkloc );
        cond21b = ( dis1<sep1 && ~cond21a );
        cond21c = ( dis1>sep3 );
        if(cond21a || cond21b || cond21c)
            if(cond21c)
                rmind = [];
                isind = ii-1;
                itvl = [Bp2.selpkloc(ii-1)+ceil(sep1),Bp2.selpkloc(ii)];
                if(ii-2>10)
                    rge = (max(ii-1-10,1):ii-2);
                else
                    rge = (1:10);
                end
            elseif(cond21b)
                adind = [];
                itvl = [Bp2.selpkloc(ii-2)+ceil(sep1),Bp2.selpkloc(ii+1)];
                if(ii-2>10)
                    rge = (max(ii-1-10,1):ii-2);
                else
                    rge = (1:10);
                end
            elseif(cond21a)
                [~,ind] = min([sum(abs(diff(Bp2.selpkloc(ii-2:ii)))),sum(abs(diff(Bp2.selpkloc(ii-1:ii+1))))]);
                if(ind==1)
                    rmind = ii;
                    isind = ii;
                    itvl = [Bp2.selpkloc(ii-1)+ceil(sep1),Bp2.selpkloc(ii+1)];
                    if(ii-1>10)
                        rge = (max(ii-10,1):ii-1);
                    else
                        rge = (1:10);
                    end
                elseif(ind==2)
                    rmind = ii-1;
                    isind = ii-1;
                    itvl = [Bp2.selpkloc(ii-2)+ceil(sep1),Bp2.selpkloc(ii-1)];
                    if(ii-2>10)
                        rge = (max(ii-1-10,1):ii-2);
                    else
                        rge = (1:10);
                    end
                end
            end
            coritvl = zeros(numel(rge),2);
            coritvlpkloc = zeros(1,numel(rge));
            for jj=1:numel(rge)
                coritvl(jj,:) = [max(Bp2.selpkloc(rge(jj))-ceil(sep0),1),min(Bp2.selpkloc(rge(jj))+floor(sep0),numel(Bp2.data))];
                coritvlpkloc(jj) = Bp2.selpkloc(rge(jj));
            end
            [pkind,sel,selflg,Bp2] = selpk2(Bp2,itvl,coritvl,coritvlpkloc,[mean(Bp2.data)+2*std(Bp2.data),mean(max([Bp2.lfslp;Bp2.rgslp],[],1))+2*std(max([Bp2.lfslp;Bp2.rgslp],[],1))],[1 1 1 1]);
            if(numel(pkind)>0)
                if(cond21c)
                    adind = pkind(find(sel==1));
                elseif(cond21b)
                    pkrmind = Bp2.pkloc(pkind(find(sel==0)));
                    for jj=[ii-1,ii]
                        if(numel(find(pkrmind==Bp2.selpkloc(ii-1)))>0)
                            rmind = ii-1;
                            isind = rmind;
                        elseif(numel(find(pkrmind==Bp2.selpkloc(ii)))>0)
                            rmind = ii;
                            isind = rmind;
                        else
                            rmind = [];
                            isind = ii;
                        end
                    end
                elseif(cond21a)
                    adind = pkind(find(sel==1));
                end
                Bp2.selpkloc = [Bp2.selpkloc(1:isind),Bp2.pkloc(adind),Bp2.selpkloc(isind+1:end)];
                Bp2.selpkloc(rmind) = [];
                ii = ii - numel(rmind) + numel(adind);
                if(cond21c==1)
                    ii = ii - 1;
                end
            end
        end
    else
        flag = 0;
    end
end

% find corresponding peak in ic
pkloc = zeros(size(Bp2.selpkloc));
for ii=1:numel(Bp2.selpkloc)
    wdt = ceil(srate/(acqfrq-0.5));
    [~,ind] = findpeaks(Fn(max(Bp2.selpkloc(ii)-wdt,1):min(Bp2.selpkloc(ii)+wdt,size(Fn,2))));
    if(numel(ind)>0)
        ind=ind+max(Bp2.selpkloc(ii)-wdt,1)-1;
        [~,ind1] = max(Fn(ind));
        pkloc(ii) = ind(ind1);
    else
        pkloc(ii) = Bp2.selpkloc(ii);
    end
end

fprintf('Multiple scale cardiac peak detection done -- %f second\n',toc(tmbp));
end



function [pkind,sel,selflg,Bp2] = selpk2(Bp2,itvl,coritvl,coritvlpkloc,thr,cond)
pkind = find(((Bp2.pkloc>=itvl(1) & Bp2.pkloc<itvl(2)))==1);
if(numel(pkind)==0)
    sel = [];
    return;
end
cond(3) = (cond(3) && numel(coritvlpkloc)>0);
pkloc = Bp2.pkloc(pkind);
selflg = zeros(numel(pkloc),numel(cond));
if( cond(1)>0 )
    [~,Ind] = max(max([Bp2.lfslp(pkind);Bp2.rgslp(pkind)],[],1));
    selflg(Ind,1) = 1;
end
if( cond(2)>0 )
    [~,Ind] = max(Bp2.pkamp(pkind));
    selflg(Ind,2) = 1;
end
if( cond(3)>0 )
    maxcor = zeros(size(coritvlpkloc));
    for jj=1:size(coritvl,1)
        [Corr,lag] = xcorr(Bp2.data(coritvl(jj,1):coritvl(jj,2)),Bp2.data(itvl(1):itvl(2)));
        [~,ind] = max(Corr);
        maxcor(jj) = itvl(1)+coritvlpkloc(jj)-coritvl(jj,1)-lag(ind);
    end
    selsc = abs(pkloc-mean(maxcor));
    [~,Ind] = min(selsc);
    selflg(Ind,3) = 1;
end
if( cond(4)>0 )
    for ii=1:numel(pkloc)
        selflg(ii,4) = ( Bp2.pkamp(pkind(ii))>=thr(1) || max(Bp2.lfslp(pkind(ii)),Bp2.rgslp(pkind(ii)))>=thr(2) );
    end
else
    selflg(:,4) = 1;
end

sel = ones(numel(pkloc),1);
for ii=1:numel(cond)
    if(cond(ii)>0)
        sel = ( sel & selflg(:,ii) );
    end
end
if(numel(find(sel==1))==0)
    if(cond(3)==1)
        sel = ( selflg(:,1) & selflg(:,2) & selflg(:,3) );
        if(numel(find(sel==1))==0)
            sel = ( selflg(:,4) & selflg(:,3) );
            if(numel(find(sel==1))==0)
                ind1 = [find(selflg(:,2)&selflg(:,3)==1),find(selflg(:,1)&selflg(:,3)==1)];
                if(numel(ind1)==0)
                    ind1 = [find(selflg(:,2)==1),find(selflg(:,1)==1)];
                end
                [~,ind] = min(selsc(ind1));
                sel(ind1(ind)) = 1;
            end
        end
    elseif(cond(3)==0)
        sel = ( selflg(:,1) & selflg(:,2) & selflg(:,4) );
        if(numel(find(sel==1))==0)
            ind1 = find(selflg(:,2)&selflg(:,4)==1);
            if(numel(ind1)==0)
                ind1 = find(selflg(:,1)&selflg(:,4)==1);
                if(numel(ind1)==0)
                    ind1 = find(selflg(:,2)==1);
                end
            end
            sel(ind1) = 1;
        end
    end
end
end



function pow = power1(Fn,srt,tmwd,tmres,bndr)
fprintf('Calculating alpha power ...\n');
tm = tic;

wsz = 512;  ovlp = 256;  nfft = 1024;
PSD = zeros(size(Fn,1),nfft/2+1);
frq = zeros(1,nfft/2+1);
for ii=1:size(Fn,1)
    [PSD(ii,:),frq(:)] = pwelch(Fn(ii,:),wsz,ovlp,nfft,srt);
end
mpsd = mean(PSD,1);
[~,loc] = max(mpsd(frq>=bndr(1)&frq<=bndr(2)));
loc = loc + numel(frq(frq<bndr(1)));
pkfq = frq(loc);
pkrg = [pkfq-2,pkfq+2];

perlg = tmwd*srt;
stlg = tmres*srt;
sttt = floor((size(Fn,2)-perlg)/stlg);
pow = zeros(size(Fn,2)/(tmres*srt),size(Fn,1));
for ii = 1:sttt+1
    spt = zeros(perlg,size(Fn,1));
    frq = (1:1:perlg)/tmwd;
    for jj = 1:size(Fn,1)
        spt(:,jj) = abs(fft(transpose(Fn(jj,stlg*(ii-1)+1:stlg*(ii-1)+perlg).*transpose(hann(perlg))))/srt).^2;
    end
    pow(ii,:) = mean(spt(frq>=pkrg(1) & frq<=pkrg(2),:),1);
end
pow(sttt+2:end,:) = repmat(pow(sttt+1,:),size(pow(sttt+2:end,:),1),1);

fprintf('EEG alpha power calculation done -- ');
toc(tm);
fprintf('\n');
end



% add new markers at the end of EEG.event
function NewEvent=AddMarker(event, MarkerTiming,duration,channel,type,code)

CurrEventSize=max(size(event,2),size(event,1));

for i=1:size(MarkerTiming,2)
    loc=CurrEventSize+i;
    event(loc).latency=MarkerTiming(i);
    event(loc).duration=duration;
    event(loc).channel=channel;
    event(loc).bvtime=[];
    event(loc).bvmknum=CurrEventSize+i;
    event(loc).type=type;
    event(loc).code=code;
    event(loc).urevent=CurrEventSize+i;
end
NewEvent=event;
end



% sort markers based on the duration timing 
function sorted = sortlatency(event,clmn2srt)
EventTable=struct2table(event);
SrtEventTable=sortrows(EventTable,clmn2srt);
NumEvent=1:1:size(SrtEventTable,1);
SrtEventTable.bvmknum=NumEvent';
SrtEventTable.urevent=NumEvent';
sorted=table2struct(SrtEventTable);
end



% identify ICs related to artifacts (updated by Kaylee and Ahmad)
function [cbicind,saccade,blink,topomap,spectrumy,tpblink,tpsac,smolregion,singchan,muscleloc] = icid(ic,A,mixsig,srte,mriperd,flg)
tol = 1.0e-6;
szic = size(ic);
szmixmat = size(A);
icnum = szic(1);
szsm = szic(2);
nbch = szmixmat(1);
[nu,de] = butter(3,0.0008,'high');
for ii=1:icnum
    ic(ii,1:szsm) = filtfilt(nu,de,ic(ii,1:szsm));
end
for ii=1:nbch
    mixsig(ii,1:szsm) = filtfilt(nu,de,mixsig(ii,1:szsm));
end

% ic statistics
[icave,icstd,~,~,ickurt] = stat(ic);
regthrs = [0.3,0.1,0.5];
grdn = 100;
thrsnum = size(regthrs,2);
diffr(1) = linintep([5,15],[1/3,0.45],mean(abs(ickurt(1:icnum)))); %Equation S10
diffr(2) = 0.25;
ignsignif = zeros(1,nbch);
ignsignif([17:19,21:24]) = 1;

% spectral analysis
mofrqind = zeros(1,icnum); % motion frequency (CB+RM motions) index
sptcond = zeros(1,icnum);
[spt,frq] = pwelchcal(ic,srte,512,256); %spt is the power spectral density curve, frq is the freq. of the PSD curve
sptsz = size(spt,2);
rmfrqrge = floor([0.5,4.5]/(srte/2)*(sptsz-1))+1; % rapid head motion freq. range
cbfrqrge = floor([2,7]/(srte/2)*(sptsz-1))+1; % cardioballistic freq. range
mofrqrge = floor([0.5,7]/(srte/2)*(sptsz-1))+1; % motion freq. range
nrfrqrge = floor([8,12]/(srte/2)*(sptsz-1))+1; % neuronal freq. range
bkfrqrge = floor([0.5,3]/(srte/2)*(sptsz-1))+1; % blink freq. range
sptscfrqrge = [1,floor(4/(srte/2)*(sptsz-1))+1]; % amplitude scale(S0) freq. range

deltaband = floor([0.5,4]/(srte/2)*(sptsz-1))+1;
thetaband = floor([4,8]/(srte/2)*(sptsz-1))+1;
alphaband = floor([8,12]/(srte/2)*(sptsz-1))+1;
betaband = floor([12,30]/(srte/2)*(sptsz-1))+1;
gammaband = floor([30,60]/(srte/2)*(sptsz-1))+1;

% easy one:
alphababy = zeros(1,icnum);
gammaflg = zeros(1,icnum);
notart = zeros(1,icnum);

% which has the highest PSD/range ratio

for kh = 1:icnum
    [wholefrqpk,wholefrqpkloc] = findpeaks(spt(kh,:));
    [maxpeak,maxpeakloc]=max(wholefrqpk);
    MaxPkLoc=wholefrqpkloc(maxpeakloc);
    
    deltaratio = sum(spt(kh,deltaband(1):deltaband(2)))/(deltaband(2)-deltaband(1));
    thetaratio = sum(spt(kh,thetaband(1):thetaband(2)))/(thetaband(2)-thetaband(1));
    alpharatio = sum(spt(kh,alphaband(1):alphaband(2)))/(alphaband(2)-alphaband(1));
    betaratio = sum(spt(kh,betaband(1):betaband(2)))/(betaband(2)-betaband(1));
    gammaratio = sum(spt(kh,gammaband(1):gammaband(2)))/(gammaband(2)-gammaband(1));
    if (alpharatio > max(max(deltaratio,thetaratio),betaratio)) || ((alphaband(1) < MaxPkLoc)&& (alphaband(2) > MaxPkLoc))
        alphababy(1,kh) = 1;
    else
        alphababy(1,kh) = 0;
    end
    
    if (gammaratio > max(max(max(deltaratio,thetaratio),betaratio),alpharatio))
        gammaflg(1,kh) = 1; %mark this one as muscle artifact
    end
    
end


sptampsc = zeros(1,icnum);
mofrqpknm = zeros(1,icnum);
mofrqpkrs = zeros(1,icnum);
cbfrqpkrs = zeros(1,icnum);
maxcbfrqpk = zeros(1,icnum);
avepowrm = zeros(1,icnum);
avepowcb = zeros(1,icnum);
rmconvwdt = zeros(1,icnum);
rmconvdpt = zeros(1,icnum);
mrfrqpkrs = zeros(1,icnum);
lgmopk = zeros(1,icnum);
condbkspt = zeros(1,icnum);
condscspt = zeros(1,icnum);
condlgmopkspt = zeros(1,icnum);
condrspspt = zeros(1,icnum);
if flg==1
    %figure; plot(frq,spt(1:31,:))
    %figure; plot(frq,spt(1:12,:))
    %figure;
end

% spectral peaks identification
% S0: sptampsc
% rnr: mrfrqpkrs
% pkrse: rcb the peak rise in the CB frequency range
for ii = 1:icnum
    [mofrqpk,mofrqpkloc] = findpeaks(spt(ii,mofrqrge(1):mofrqrge(2))); % find the spectrum peaks in MO freq. range in the PSD curve
    mofrqpknm(ii) = numel(mofrqpk); % how many MO peaks are there in each IC
    if(mofrqpknm(ii)>0) % If there is a peak...
        mofrqpkloc = mofrqpkloc + mofrqrge(1) - 1; % The MO peak loc + lower bound of MO freq. range - 1 (correct the peak location)
    end
    [nrfrqpk,nrfrqpkloc] = max(spt(ii,nrfrqrge(1):nrfrqrge(2))); % find the spectrum max in neuronal freq. range in the PSD curve
    nrfrqpkloc = nrfrqpkloc + nrfrqrge(1) - 1; % The NR peak loc + lower bound of NR freq. range - 1 (correct the peak location)
    [pkloc,~,~,~,dpt,~,~] = fninfo(spt(ii,1:sptsz)); % 0.05*S0
    sptampsc(ii) = max(spt(ii,sptscfrqrge(1):sptscfrqrge(2)))-min(spt(ii,sptscfrqrge(1):sptscfrqrge(2))); % S0
    minpkamp = 0.05*sptampsc(ii); % 0.05*S0
    dscdpk1 = zeros(1,mofrqpknm(ii)); % discard peak locations; zeros the size of number of peaks found in the MO freq. range
    for jj=1:mofrqpknm(ii) %1 to the number of peaks
        dpt1 = dpt(pkloc==mofrqpkloc(jj)); %where Vr is below 8Hz
        if(dpt1<minpkamp) %if the value is less than 0.05*S0...
            dscdpk1(jj) = 1; %mark the location to be discarded
        end
    end
    mofrqpknm(ii) = mofrqpknm(ii) - sum(dscdpk1); %update the number of MO peaks after some have been discarded
    mofrqpkloc = mofrqpkloc(dscdpk1==0); %update the locations of MO peaks after some have been discarded
    
    [~,minimaloc] = findpeaks(-spt(ii,1:sptsz)); %find location of all local minima of the PSD curve
    [nrpklftmin,~] = min(abs(minimaloc(minimaloc < nrfrqrge(1))-nrfrqrge(1))); % The difference between the start of NR range and the last minimal before starting the NR range
    %Table S1 aEREMCOR Peak rise in NR range
    if(numel(nrpklftmin)>0) %if the NR peak minimum exists then..
        nrpklftmin = nrfrqrge(1) - nrpklftmin; %loc of the left minima from NR
        mrfrqpkrs(ii) = nrfrqpk - min(spt(ii,nrpklftmin:nrfrqpkloc)); %if Vl exists
    elseif(numel(nrpklftmin)==0)%if Vl doesn't exist
        mrfrqpkrs(ii) = nrfrqpk - min(spt(ii,nrfrqrge(1):nrfrqpkloc));
    end
    maxcbfrqpk(ii) = min(spt(ii,cbfrqrge(1):cbfrqrge(2))); %minimum value found in the CB freq. range (in cb range, location needs to be corrected since it is before shifting)
    nbkpkrs = 0;
    nbkpkflg = 0;
    cbpkflg = 0;
    rgtminflg = 0;
    for jj=1:mofrqpknm(ii) %for 1:number of MO peaks
        [rgtmin,~] = min(abs(minimaloc(minimaloc > mofrqpkloc(jj)) - mofrqpkloc(jj))); %right min; smallest distance above the local minima from the MO peak
        [lftmin,~] = min(abs(minimaloc(minimaloc < mofrqpkloc(jj)) - mofrqpkloc(jj))); %left min; smallest distance below the local minima from the MO peak
        if( numel(lftmin)==0 ) %if no left min
            lftmin = 1; %make left min = 1
            lftminflg = 0;
        else
            lftmin = mofrqpkloc(jj) - lftmin; %Adjust the correct loc. left min is the loc, after removing left min from the MO peak
            lftminflg = 1;
        end
        rgtmin = mofrqpkloc(jj) + rgtmin; %Adjust the correct loc.right min is the loc after adding right min to the MO peak
        if( rgtminflg==0 && frq(rgtmin)<=7 ) %if local minimum is below 7 Hz
            rgtminflg = 1;
        else
            rgtminflg = 0;
        end
        
        lftrse = spt(ii,mofrqpkloc(jj))-spt(ii,lftmin); %left rise; PSD value of MO peak - left local minima
        rgtrse = spt(ii,mofrqpkloc(jj))-spt(ii,rgtmin); %right rise; PSD value of MO peak - right local minima
        pkrse = 0.5*(lftrse+rgtrse); %peak rise
        if( lftminflg==1 && frq(mofrqpkloc(jj))>0.5 && frq(mofrqpkloc(jj))<4.5 ) % If the motion peak is in the RM range
            if( min(lftrse,rgtrse)>minpkamp )
                nbkpkrs = max( nbkpkrs, pkrse );
                nbkpkflg = 1;
            end
        end
        if( lftminflg==1 && frq(mofrqpkloc(jj))>2 && frq(mofrqpkloc(jj))<7 ) % If the motion peak is in the CB range
            if( pkrse>0.2*sptampsc(ii) ) % rcb>0.2*S0 condition (3) Large MO peak
                cbfrqpkrs(ii) = max(cbfrqpkrs(ii),pkrse);% max(rcb)>0.2S0
                maxcbfrqpk(ii) = max(maxcbfrqpk(ii),spt(ii,mofrqpkloc(jj))); %????
                cbpkflg = 1;
            end
        end
        mofrqpkrs(ii) = max(mofrqpkrs(ii),pkrse); %max peak rise in cb range
    end
    avepowrm(ii) = mean(spt(ii,rmfrqrge(1):rmfrqrge(2)))-min(spt(ii,1:nrfrqpkloc)); % mean(Srm)-Smin
    avepowcb(ii) = mean(spt(ii,cbfrqrge(1):cbfrqrge(2)))-min(spt(ii,1:nrfrqpkloc)); % mean(Scb)-Smin
    
    % spectrum slope in motion freq range
    rmsptsc = max(spt(ii,rmfrqrge(1):rmfrqrge(2)))-min(spt(ii,rmfrqrge(1):rmfrqrge(2)));
    rmsptnrm = ( spt(ii,1:sptsz) - min(spt(ii,rmfrqrge(1):rmfrqrge(2))) ) / rmsptsc; % Normalize the spectrum in RM freq. range
    deri1 = diff(rmsptnrm(1:sptsz)); % the first derivative of normalized spectrum
    smslp = deri1;
    % smooth the local fluctuation of each point by taking the average
    % with its neighbors
    for jj=2:sptsz-2
        smslp(jj) = mean(deri1(jj-1:jj+1));
    end
    deri2 = diff(smslp); % the second derivative of normalized spectrum
    trnptflg = zeros(size(deri2));
    for jj=2:sptsz-3
        if(deri2(jj-1)>=0 && deri2(jj+1)<=0 && deri2(jj-1)-deri2(jj+1)>0.02) % S"(v0-)>=0, S"(V0+)<=0 and S"(V0-)-S"(V0+)>0.02
            trnptflg(jj) = 1;
        end
    end
    trnpt = find(trnptflg==1); % Find reflection points
    % Find reflection points in RM range
    trnpt = trnpt(trnpt>rmfrqrge(1));
    trnpt = trnpt(trnpt<rmfrqrge(2));
    for jj=1:numel(trnpt)
        if(jj<numel(trnpt))
            rgtpt = trnpt(jj+1);
        elseif(jj==numel(trnpt)) %last reflection point
            [~,rgtpt] = min(rmsptnrm(trnpt(jj):rmfrqrge(2)));
            rgtpt = rgtpt + trnpt(jj) - 1;
        end
        tmpwdt = (rgtpt-trnpt(jj))*(srte/2)/(sptsz-1);
        tmpdpt = rmsptnrm(trnpt(jj))-rmsptnrm(rgtpt); %Rrm'=S(V0)-S(V'0)
        if( tmpdpt>rmconvdpt(ii) )
            rmconvwdt(ii) = tmpwdt;
            rmconvdpt(ii) = tmpdpt;
        end
    end
    
    % spectrum slope in blink freq range
    bksptsc = max(spt(ii,bkfrqrge(1):bkfrqrge(2))) - min(spt(ii,bkfrqrge(1):bkfrqrge(2)));
    bksptnrm = ( spt(ii,bkfrqrge(1):bkfrqrge(2)) - min(spt(ii,bkfrqrge(1):bkfrqrge(2))) ) / bksptsc;
    bkflt = std(diff(diff(bksptnrm)));
    nbkpkind = 0;
    if(nbkpkflg==1 && nbkpkrs>2/3*sptampsc(ii))
        nbkpkind = 1;
    end
    condbkspt(ii) = or( bkflt<0.2 && not(nbkpkind), mofrqpknm(ii)==0 );
    condscspt(ii) = condbkspt(ii);
    
    % selection
    sussptcond = 0;
    if(mofrqpknm(ii)>=1)
        sussptcond = ( lftminflg==0 && rgtminflg==0 );
    end
    minrmfrqval = min(spt(ii,rmfrqrge(1):rmfrqrge(2)));
    avemopkval = mean(spt(ii,mofrqpkloc));
    lgmopk(ii) = ( mofrqpkrs(ii)>0.2*sptampsc(ii) );
    smmopk = ( not(lgmopk(ii)) && avemopkval>nrfrqpk-2 ); %Condition (4) Small MO peaks
    smnrpk = or( mrfrqpkrs(ii)<0.2*sptampsc(ii), minrmfrqval>nrfrqpk-2 ); %Condition (5) Small NR peak
    domnrpk = ( mrfrqpkrs(ii)>sptampsc(ii)/3 ); % If rise in NR range is larger than 0.33 of the largest rise in all spectrum range
    if(domnrpk==1)
        cond1 = ( cbpkflg==1 && cbfrqpkrs(ii)>mrfrqpkrs(ii)-3 ); %Condition (iii):max(rcb)>rnr-3
        cond2 = ( cbpkflg==1 && maxcbfrqpk(ii)>nrfrqpk-3 && avepowcb(ii)>mrfrqpkrs(ii)/3 ); % Condition (iv)maxcbfrqpk: max(Scb) and nrfrqpk: max(Snr)
        cond3 = ( avepowrm(ii)>1.75*mrfrqpkrs(ii) ); %condition (ii) in aEremcor mean(S(v)-Smin>a2rnr
        accnrpk = ( cond1 || cond2 || cond3 );
    else
        accnrpk = 1;
    end

    if( accnrpk && not(sussptcond) )
        if( mofrqpknm(ii)>=1 && lgmopk(ii) )
            mofrqind(ii) = 1;
            if( cbpkflg==1 )
                condlgmopkspt(ii) = 1;
                %               if (alphababy(ii) == 1)
                %                   condlgmopkspt(ii) = 0;
                %               end
            end
        elseif( mofrqpknm(ii)>=1 && and(smmopk,smnrpk) )
            mofrqind(ii) = 1;
        elseif( rmconvdpt(ii)>0.15 && smnrpk )  %Condition 6: obvious negative convexity
            mofrqind(ii) = 1;
        end
    end
    condrspspt(ii) = accnrpk;
end

% topo map analysis
motpind = zeros(1,icnum);
polcond = zeros(1,icnum);
cbtpind = zeros(1,icnum);
Bipflg = zeros(1,icnum); %Bipolar right and left flag
[xq, yq] = meshgrid(-1:2/grdn:1, -1:2/grdn:1);
regnm = zeros(icnum,thrsnum,2);
regArcnm = zeros(icnum,thrsnum,2);
regA = zeros(icnum,thrsnum,2);
regArcA = zeros(icnum,thrsnum,2);
regcntval = zeros(icnum,thrsnum,2);
regcnt2org = zeros(icnum,thrsnum);
reg12cnt2org = zeros(icnum,1);
regX = zeros(icnum,thrsnum,2);
regY = zeros(icnum,thrsnum,2);
regR = zeros(icnum,thrsnum,2);
regT = zeros(icnum,thrsnum,2);
rregnm = zeros(icnum,thrsnum);
rregArcnm = zeros(icnum,thrsnum);
rregA = zeros(icnum,thrsnum);
rregArcA = zeros(icnum,thrsnum);
condbktp = zeros(1,icnum);
condsctp = zeros(1,icnum);
condrmmot = zeros(1,icnum);
sustpcond = zeros(1,icnum);
x = zeros(1,nbch);
y = zeros(size(x));
chpos = [[-90,-72];[90,72];[-60,-51];[60,51];[-45,0];[45,0];[-60,51];[60,-51];...
    [-90,72];[90,-72];[-90,-36];[90,36];[-90,0];[90,0];[-90,36];[90,-36];...
    [45,90];[0,0];[45,-90];[90,-90];[-31,-46];[31,46];[-31,46];[31,-46];...
    [-69,-21];[69,21];[-69,21];[69,-21];[-113,18];[113,-18];[67,-90]];

for ii=1:nbch
    x(ii) = chpos(ii,1)*pi/180 * cos(chpos(ii,2)*pi/180);
    y(ii) = chpos(ii,1)*pi/180 * sin(chpos(ii,2)*pi/180);
end
x1 = x./max(sqrt(x.^2+y.^2));
y1 = y./max(sqrt(x.^2+y.^2));
x = x1;
y = y1;
toosmall = zeros(2,size(ic,1)); %maybe make this size (ii,jj,kk)?
pertop_top = zeros(size(ic,1),3);
pertop_bot = zeros(size(ic,1),3);
%   perbot = zeros(size(ic,1),3);
for ii=1:icnum
    v(1:nbch) = A(1:nbch,ii) / max(abs(A(1:nbch,ii)));
    vq = griddata(x,y,v,xq,yq,'v4');
    vq = vq / max(max(abs(vq)));
    vq(xq.^2 + yq.^2 > 1) = 0; %Topographic map in a circle imagesc(vq)
    vqarc = vq;
    %if flg==1
    %    subplot(8,4,ii); imagesc(vq)
    %end
    %imagesc(vq)
    vqarc(sqrt(xq.^2 + yq.^2) < 0.8) = 0; %Map boundry of 0.2 (arc ciecle)
    if(ii==1)
        bdmap = zeros(size(vq));
        bdinfo = bwboundaries(vq,8,'noholes');
        for ll=1:length(bdinfo{1})
            bdmap(bdinfo{1}(ll,1),bdinfo{1}(ll,2)) = 1; %bdmap is the outer circle of the map
        end
    end
    totA = numel(vq(xq.^2+yq.^2<=1)); %The number of pixels in outer circle
    totArcA = numel(vqarc(vqarc~=0)); %The number of pixels in the arc (the area between outer circle and inner circle)
    
    % regions definition
    for jj=1:thrsnum % jj = 1,2 are thresholds for the polarity regions (page 3 of supplementary for aE-REMOCOR); jj=1 (K=0.3) are primary polarity regions, jj=2 (K'=0.1) are secondary polarity regions
        map = ones(size(vq));

        squaretopmap = zeros(size(vq));
        squarebotmap = zeros(size(vq));
        emptymap = ones(size(vq));
        emptymap(xq.^2 + yq.^2 > 1) = 0;
        squaretopmap(1:20,20:80)=1;
        squarebotmap(80:100,20:80)=1;
        mapbaby_top = emptymap.*squaretopmap;
        mapbaby_bot = emptymap.*squarebotmap;
        mapbaby = mapbaby_top + mapbaby_bot;
        
        squarerightmap = zeros(size(vq));
        squareleftmap = zeros(size(vq));
        squarerightmap(20:80,80:100)=1;
        squareleftmap(20:80,0.01:22)=1;
        mapbaby_right = emptymap.*squarerightmap;
        mapbaby_left = emptymap.*squareleftmap;
        mapbaby_bip = mapbaby_right + mapbaby_left;
        %figure;
        %imagesc(mapbaby_bip)

        
        squarerightangmap = zeros(size(vq));
        squareleftangmap = zeros(size(vq));
        squarerightangmap(0.01:30,0.01:30)=1;
        squareleftangmap(70:100,70:100)=1;
        mapbaby_rightang = emptymap.*squarerightangmap;
        mapbaby_leftang = emptymap.*squareleftangmap;
        mapbaby_ang = mapbaby_rightang + mapbaby_leftang;
        
        
        
        squarerightangmap1 = zeros(size(vq));
        squareleftangmap1 = zeros(size(vq));
        squarerightangmap1(70:100,0.01:30)=1;
        squareleftangmap1(0.01:30,70:100)=1;
        mapbaby_rightang1 = emptymap.*squarerightangmap1;
        mapbaby_leftang1 = emptymap.*squareleftangmap1;
        mapbaby_ang1 = (mapbaby_rightang1 + mapbaby_leftang1);
        
        
        squarerightangmap2 = zeros(size(vq));
        squareleftangmap2 = zeros(size(vq));
        squarerightangmap2(0.01:30,0.01:30)=1;
        squareleftangmap2(0.01:30,70:100)=1;
        mapbaby_rightang2 = emptymap.*squarerightangmap2;
        mapbaby_leftang2 = emptymap.*squareleftangmap2;
        mapbaby_ang2 = (mapbaby_rightang2 + mapbaby_leftang2);
        
        
        
        
        
        map(xq.^2 + yq.^2 > 1) = 0;
        map(vq>regthrs(jj))=0; %positive polarity
        map(vq<-regthrs(jj))=0; %negative polarity
        maplabl = bwlabel(map,8); %separate founded regions (neutral ones |g(r)|<thrdnum(K)
        selregR = sign(maplabl);
        rregnm(ii,jj) = max(max(maplabl)); % Count number of brain regions for each IC(ii) and each thrsnum(jj); column 1 is # of primary polarity regions, column 2 is # of secondary polarity regions, column 3 is # of other regions?
        npatA = zeros(1,rregnm(ii,jj));
        dscdR = 0;
        ignR = 0;
        for ll=1:rregnm(ii,jj)
            npatA(ll) = sum(sign(maplabl(maplabl==ll)))/totA; %The percentage of neutral area to the area of outer circle
            if(npatA(ll)<=0.04 && npatA(ll)>0.001 && dscdR<2) % If the region is very small, just discard it
                dscdR = dscdR + 1; % Add it to Discarded neutral  brain regions
                selregR(maplabl==ll) = 0;
            end
            if(npatA(ll)<=0.001) % If the region is very small, just discard it
                ignR = ignR + 1;
                selregR(maplabl==ll) = 0;
            end
        end
        rregnm(ii,jj) = rregnm(ii,jj) - dscdR - ignR; %number of neutral regions after discarding/ignoring
        map(selregR==0) = 0; %removes the discard/ignore neutral regions from the map
        arcmap = map;
        arcmap(vqarc==0) = 0; %overlap between neutral region map and the arc
        rregArcnm(ii,jj) = max(max(bwlabel(arcmap,8)));
        rregA(ii,jj) = numel(map(map~=0))/totA;  %The percentage of neutral area to the area of outer circle
        rregArcA(ii,jj) = numel(arcmap(arcmap~=0))/totArcA;  %The percentage of neutral area in the arc to the area of outer circle
        for kk=1:2 %1 is positive and 2 negative regions
            ignmap = zeros(size(vq));
            ignmap((-1)^(kk+1)*vq<=regthrs(jj)) = 1;
            ignmap(xq.^2 + yq.^2 > 1) = 0;
            ignmaplabl = bwlabel(ignmap,8);
            ignpatA = zeros(1,max(max(ignmaplabl)));
            for ll=1:max(max(ignmaplabl))
                ignpatA(ll) = sum(sign(ignmaplabl(ignmaplabl==ll)))/totA; %percentage ratio of anything outside the ll region
                %discarding small regions
                if(ignpatA(ll)<=0.025 && ignpatA(ll)>0.001)
                    ignmaplabl(ignmaplabl==ll) = 0;
                end
                if(ignpatA(ll)<=0.001)
                    ignmaplabl(ignmaplabl==ll) = 0;
                end
            end
            map = vq;
            map((-1)^(kk+1)*vq>regthrs(jj)) = 1;
            map((-1)^(kk+1)*vq<=regthrs(jj)) = 0;
            maplabl = bwlabel(map,8);
            selreg = sign(maplabl);  %selected region (one of the positive or negative region)
            
            together_top=mapbaby_top.*selreg;
            together_bot=mapbaby_bot.*selreg;
            pertop_top(ii,jj,kk) = sum(sum(together_top))/sum(sum(mapbaby_top));
            pertop_bot(ii,jj,kk) = sum(sum(together_bot))/sum(sum(mapbaby_bot));
            together_bip=mapbaby_bip.*selreg;
            pertop_bip(ii,jj,kk) = sum(sum(together_bip))/sum(sum(mapbaby_bip));
            
            together_ang=mapbaby_ang.*selreg;
            pertop_ang(ii,jj,kk) = sum(sum(together_ang))/sum(sum(mapbaby_ang));
            
            together_ang1=mapbaby_ang1.*selreg;
            pertop_ang1(ii,jj,kk) = sum(sum(together_ang1))/sum(sum(mapbaby_ang1));
            
            together_ang2=mapbaby_ang2.*selreg;
            pertop_ang2(ii,jj,kk) = sum(sum(together_ang2))/sum(sum(mapbaby_ang2));
            % else
            % together=mapbaby.*selreg;
            %pertop(ii,jj,kk) = sum(sum(together))/sum(sum(mapbaby));
            %end
            
            regnm(ii,jj,kk) = max(max(maplabl));
            patA = zeros(1,regnm(ii,jj,kk));
            dscd = 0;
            ign = 0;
            arcmap = sign(abs(vqarc)); %part of this region that is inside the arc
            arcmap(ignmaplabl>0) = 0;
            
            for ll=1:regnm(ii,jj,kk)
                patA(ll) = sum(sign(maplabl(maplabl==ll)))/totA; % The percentage ratio of the ll region over the area f outer circle
                %If it's too small, ignore it
                if(patA(ll)<=0.025 && patA(ll)>0.001 && dscd<2)
                    if kk==1;
                        toosmall (1,ii) = toosmall(1,ii)+patA(ll);
                    else
                        toosmall (2,ii) = toosmall(2,ii)+patA(ll);
                    end
                    dscd = dscd + 1;
                    selreg(maplabl==ll) = 0;
                    arcmap(maplabl==ll) = 0;
                end
                if(patA(ll)<=0.001)
                    ign = ign + 1;
                    selreg(maplabl==ll) = 0;
                    arcmap(maplabl==ll) = 0;
                end
            end
            regnm(ii,jj,kk) = regnm(ii,jj,kk) - dscd - ign; %number of neutral regions after discarding/ignoring
            arcmap(selreg==0) = 0;
            regArcnm(ii,jj,kk) = max(max(bwlabel(arcmap,8)));
            
            regA(ii,jj,kk) = numel(vq(selreg>0))/totA; %percentage of the positive/negative region area to the entire circle area
            regArcA(ii,jj,kk) = numel(vqarc((-1)^(kk+1)*vqarc>regthrs(jj)&selreg>0))/totArcA; %percentage of positive/negative region arc area to the entire arc area
            regX(ii,jj,kk) = mean(xq(selreg>0));
            regY(ii,jj,kk) = mean(yq(selreg>0));
            regR(ii,jj,kk) = mean(sqrt((xq(selreg>0)).^2+(yq(selreg>0)).^2)); % The radius of the center of positive/negative region
            regT(ii,jj,kk) = atan(abs(regY(ii,jj,kk))/abs(regX(ii,jj,kk)));
            % Fixing the Theta
            if(regX(ii,jj,kk)<0 && regY(ii,jj,kk)>0)
                regT(ii,jj,kk) = pi - regT(ii,jj,kk);
            elseif(regX(ii,jj,kk)<0 && regY(ii,jj,kk)<0)
                regT(ii,jj,kk) = regT(ii,jj,kk) + pi;
            elseif(regX(ii,jj,kk)>0 && regY(ii,jj,kk)<0)
                regT(ii,jj,kk) = 2*pi - regT(ii,jj,kk);
            end
            if(regArcA(ii,jj,kk)>tol)
                regX(ii,jj,kk) = regR(ii,jj,kk)*cos(regT(ii,jj,kk));
                regY(ii,jj,kk) = regR(ii,jj,kk)*sin(regT(ii,jj,kk));
            end
            [~, locx] = min(abs((-1:2/grdn:1)-regX(ii,jj,kk)));
            [~, locy] = min(abs((-1:2/grdn:1)-regY(ii,jj,kk)));
            regcntval(ii,jj,kk) = vq((locx-1)*(grdn+1)+locy); %find the gradient value at the x,y location
        end
        regcnt2org(ii,jj) = lne2pt([regX(ii,jj,1) regY(ii,jj,1) 0],[regX(ii,jj,2) regY(ii,jj,2) 0],[0 0 0]); %finds distance d
    end
    if flg==1
        %% blink identification
        bkthrsapp = 1;
        condbktp1=0;
        condbktp2=0;
        condbktp3=0;
        [~, nseposx] = min(abs((-1:2/grdn:1)-0));
        [~, nseposy] = min(abs((-1:2/grdn:1)-1));
        nseposval = vq((nseposx-1)*(grdn+1)+nseposy);
        for kk=1:2
            if( (-1)^(kk+1)*nseposval>2/3 )
                bkregsgn = kk;
                map = vq;
                map((-1)^(bkregsgn+1)*vq>regthrs(bkthrsapp)) = 1;
                map((-1)^(bkregsgn+1)*vq<=regthrs(bkthrsapp)) = 0;
                mapnse = bwselect(map,nseposx,nseposy,8);
                arcmapnse = mapnse;
                arcmapnse(vqarc==0) = 0;
                arcmapnsenm = max(max(bwlabel(arcmapnse,8)));
                smapnse = mapnse;
                smapnse(abs(vq)<2/3) = 0;
                smapnsenm = max(max(bwlabel(smapnse,8)));
                sArcmapnse = smapnse;
                sArcmapnse(vqarc==0) = 0;
                bkArcA = sum(sum(arcmapnse==1))/totArcA;
                bkA = sum(sum(mapnse==1))/totA;
                bksArcA = sum(sum(sArcmapnse==1))/totArcA;
                bksA = sum(sum(smapnse==1))/totA;
                bkbdinfo = bwboundaries(mapnse,8,'noholes');
                bkbd = max(length(bkbdinfo{1})/floor(grdn/2)-2*pi*bkArcA,0);
                rge = [min(((-1)^bkregsgn)*[-0.3,0.45]),max(((-1)^bkregsgn)*[-0.3,0.45])];
                nbkA = numel(vq(xq.^2+yq.^2<=1 & vq>=rge(1) & vq<=rge(2)))/totA;
                bkregX = mean(xq(mapnse==1));
                bkregY = mean(yq(mapnse==1));
                bkregR = mean(sqrt((xq(mapnse==1)).^2+(yq(mapnse==1)).^2));
                bkregT = atan(abs(bkregY)/abs(bkregX));
                if(bkregX<0 && bkregY>0)
                    bkregT = pi - bkregT;
                elseif(bkregX<0 && bkregY<0)
                    bkregT = bkregT + pi;
                elseif(bkregX>0 && bkregY<0)
                    bkregT = 2*pi - bkregT;
                end
                if(bkArcA>tol)
                    bkregX = bkregR*cos(bkregT);
                    bkregY = bkregR*sin(bkregT);
                end
                
                Anglim = pi/2 - 4/9*pi;
                [~, loc1x] = min(abs((-1:2/grdn:1)-cos(Anglim)));
                [~, loc1y] = min(abs((-1:2/grdn:1)-sin(Anglim)));
                [~, loc2x] = min(abs((-1:2/grdn:1)+cos(Anglim)));
                [~, loc2y] = min(abs((-1:2/grdn:1)-sin(Anglim)));
                condbktp1 = ( arcmapnsenm==1 && smapnsenm==1 && bkregY>0 && abs(bkregX)<bkregY*tan(35*pi/180) && smapnse(loc1y,loc1x)==0 && smapnse(loc2y,loc2x)==0 );
                if(condbktp1)
                    secAlim = zeros(1,2);
                    condbktp2a = ( bkArcA>13/72 && bkArcA<4/9 && bksArcA>0.45*bkArcA );
                    SectorArea = ( 2*pi*bkArcA-sin(2*pi*bkArcA) )/2 /pi;
                    secAlim(1:2) = [2,0.9].*( 2*pi*[13/72,4/9]-sin(2*pi*[13/72,4/9]) )/2 /pi;
                    Alim1 = max(0.8*SectorArea, secAlim(1));
                    Alim2 = min(SectorArea+0.8*sin(bkArcA*2*pi)/2/pi,secAlim(2));
                    condbktp2b = ( bkA>Alim1 && bkA<Alim2 && bksA>0.25*bkA );
                    sectbd = 1.75 * 2*sin(pi*bkArcA);
                    condbdtp2c = ( bkbd<sectbd );
                    contbdtp2d = ( nbkA>0.6 || and(nbkA<0.6,not(condlgmopkspt(ii))) );
                    condbktp2 = ( condbktp2a && condbktp2b && condbdtp2c && contbdtp2d );
                end
                [~, locx] = min(abs((-1:2/grdn:1)-regX(ii,bkthrsapp,kk)));
                [~, locy] = min(abs((-1:2/grdn:1)-regY(ii,bkthrsapp,kk)));
                condbktp3 = ( mapnse(locy,locx)==1 );
            end
        end
        if( condbktp1 && condbktp2 && condbktp3 )
            condbktp(ii) = 1;
        end
        if( rregnm(ii,3)<=1 && all(eq(regArcnm(ii,3,1:2),1)) && all(eq(regnm(ii,3,1:2),1)) )
            [~,minind] = min(regA(ii,3,1:2));
            maxind = find([1,2]~=minind);
            sustpcond(ii) = ( regY(ii,3,minind)>0 && regY(ii,3,maxind)<0 && ...
                max(regA(ii,3,1:2))>0.45 && min(regA(ii,3,1:2))<0.25 && ...
                abs(regX(ii,3,minind))<abs(regY(ii,3,minind))*tan(35*pi/180) );
        end
        
        %% saccade identification
        scthrsapp = 1;
        condsctp1=0;
        condsctp2=0;
        condsctp3=0;
        [~, scpos1x] = min(abs((-1:2/grdn:1)-(-0.9*cos(35/180*pi))));
        [~, scpos1y] = min(abs((-1:2/grdn:1)-0.9*sin(35/180*pi)));
        scpos1val = vq((scpos1x-1)*(grdn+1)+scpos1y);
        [~, scpos2x] = min(abs((-1:2/grdn:1)-0.9*cos(35/180*pi)));
        [~, scpos2y] = min(abs((-1:2/grdn:1)-0.9*sin(35/180*pi)));
        scpos2val = vq((scpos2x-1)*(grdn+1)+scpos2y);
        scArcA = zeros(1,2);
        scA = zeros(1,2);
        scsArcA = zeros(1,2);
        scsA = zeros(1,2);
        scbd = zeros(1,2);
        scregsgn = zeros(1,2);
        scregX = zeros(1,2);
        scregY = zeros(1,2);
        scregR = zeros(1,2);
        scregT = zeros(1,2);
        for kk=1:2
            if( (-1)^(kk+1)*scpos1val>regthrs(scthrsapp) && (-1)^(kk+1)*scpos2val<-regthrs(scthrsapp) )
                scregsgn(1) = kk;
                scregsgn(2) = find([1,2]~=scregsgn(1));
                map = vq;
                map((-1)^(scregsgn(1)+1)*vq>regthrs(scthrsapp)) = 1;
                map((-1)^(scregsgn(1)+1)*vq<=regthrs(scthrsapp)) = 0;
                mapsc1 = bwselect(map,scpos1x,scpos1y,8);
                map((-1)^(scregsgn(2)+1)*vq>regthrs(scthrsapp)) = 1;
                map((-1)^(scregsgn(2)+1)*vq<=regthrs(scthrsapp)) = 0;
                mapsc2 = bwselect(map,scpos2x,scpos2y,8);
                nscA = numel(vq(xq.^2+yq.^2<=1 & vq>=-0.3 & vq<=0.3))/totA;
                for kkp=1:2
                    if(kkp==1)
                        mapsc = mapsc1;
                    elseif(kkp==2)
                        mapsc = mapsc2;
                    end
                    scArcA(kkp) = numel(vqarc(vqarc~=0 & mapsc==1))/totArcA;
                    scA(kkp) = sum(sign(abs(vq(mapsc==1))))/totA;
                    scsArcA(kkp) = sum(sign(abs(vqarc(vqarc~=0 & mapsc==1 & abs(vq)>0.5))))/totArcA;
                    scsA(kkp) = sum(sign(abs(vq(mapsc==1 & abs(vq)>0.5))))/totA;
                    scbdinfo = bwboundaries(mapsc,8,'noholes');
                    scbd(kkp) = max(length(scbdinfo{1})/floor(grdn/2)-2*pi*scArcA(kkp),0);
                    if(scsA(kkp)>0)
                        scregX(kkp) = mean(xq(mapsc==1 & abs(vq)>0.5));
                        scregY(kkp) = mean(yq(mapsc==1 & abs(vq)>0.5));
                    else
                        scregX(kkp) = mean(xq(mapsc==1));
                        scregY(kkp) = mean(yq(mapsc==1));
                    end
                    scregR(kkp) = mean(sqrt((xq(mapsc==1)).^2+(yq(mapsc==1)).^2));
                    scregT(kkp) = atan(abs(scregY(kkp))/abs(scregX(kkp)));
                    if(scregX(kkp)<0 && scregY(kkp)>0)
                        scregT(kkp) = pi - scregT(kkp);
                    elseif(scregX(kkp)<0 && scregY(kkp)<0)
                        scregT(kkp) = scregT(kkp) + pi;
                    elseif(scregX(kkp)>0 && scregY(kkp)<0)
                        scregT(kkp) = 2*pi - scregT(kkp);
                    end
                    if(scArcA(kkp)>tol)
                        scregX(kkp) = scregR(kkp)*cos(scregT(kkp));
                        scregY(kkp) = scregR(kkp)*sin(scregT(kkp));
                    end
                end
                
                dis = lne2pt([scregX(1),scregY(1),0],[scregX(2),scregY(2),0],[0,0,0]);
                incnt = atan(abs(scregY(2)-scregY(1))/abs(scregX(2)-scregX(1)));
                condsctp1 = ( dis>0.3 && incnt<35/180*pi );
                if(condsctp1)
                    condsctp2a = ( min(scArcA)<0.3 && max(scArcA)<0.4 && sum(scArcA)<2/3 );
                    condsctp2b = ( max(scA)<0.4 && max(scsA./scA)>0.3 && all(gt(scsArcA./scArcA,0.9*scsA./scA)+eq(scsArcA./scArcA,0.9*scsA./scA)) );
                    sectbd = 4*sin(pi*scArcA);
                    condsctp2c = ( scbd(1)<sectbd(1) && scbd(2)<sectbd(2) );
                    condsctp2d = ( nscA>0.45 );
                    condsctp2e = ( abs(scA(1)-scA(2))/max(scA)<2/3 );
                    condsctp2 = ( condsctp2a && condsctp2b && condsctp2c && condsctp2d && condsctp2e);
                end
                [~, locx] = min(abs((-1:2/grdn:1)-regX(ii,scthrsapp,scregsgn(1))));
                [~, locy] = min(abs((-1:2/grdn:1)-regY(ii,scthrsapp,scregsgn(1))));
                condsctp3a = ( mapsc1(locy,locx)==1 );
                [~, locx] = min(abs((-1:2/grdn:1)-regX(ii,scthrsapp,scregsgn(2))));
                [~, locy] = min(abs((-1:2/grdn:1)-regY(ii,scthrsapp,scregsgn(2))));
                condsctp3b = ( mapsc2(locy,locx)==1 );
                condsctp3 = ( condsctp3a && condsctp3b );
            end
        end
        if( condsctp1 && condsctp2 && condsctp3 )
            condsctp(ii) = 1;
        end
    end
    %%
    % selection
    cond5 = 0;
    for kk=1:2
        kkp = find([1,2]~=kk);
        %if ((pertop_bot(ii,1,kk)>=0.35||(pertop_top(ii,1,kk)>=0.41)) && ((pertop_bot(ii,1,kkp)>=0.35)||(pertop_top(ii,1,kkp)>=0.41)))
        if ((pertop_top(ii,1,kk))>0.40) ||((pertop_top(ii,1,kkp))>0.40)||((pertop_top(ii,1,kk)+(pertop_top(ii,1,kkp)))>0.91)
            %if alphababy(ii)
            cond5 = 1;
            break;
            %end
        end
    end
    cond6 = 0;
    for kk=1:2
        kkp = find([1,2]~=kk);
        if (pertop_bip(ii,1,kk)>=0.40 && pertop_bip(ii,1,kkp)>=0.40)||(pertop_bip(ii,2,kk)>=0.41 && pertop_bip(ii,2,kkp)>=0.41)
            cond6 = 1;
            break;
        end
    end

    cond1 = ( rregnm(ii,1)<=1 ); %Condition 1: no more than one neutral region in the topographic map
    % cond2a: only one positive/negative arc and and region
    cond2a = ( all(eq(regArcnm(ii,1,1:2),1)) && all(eq(regnm(ii,1,1:2),1)) ); %Condition 2: Only one positive/negative polarity region and polarity arc region is allowed
    cond2b = 0;
    for kk=1:2 %for two polarity regions
        kkp = find([1,2]~=kk);
        if( regArcnm(ii,1,kk)==1 && regnm(ii,1,kk)==1 && regArcnm(ii,2,kkp)==1 && regnm(ii,2,kkp)==1 && regA(ii,1,kk)>regA(ii,1,kkp) )
            cond2b = 1;
            break;
        end
    end
    %     cond2 = ( cond2a || and(cond2b,condlgmopkspt(ii)) );
    [~,ind1] = max(regA(ii,1,1:2)); %find the larger region
    [~,ind2] = min(regA(ii,1,1:2)); %find the smaller region
    cond2a = ( or(regnm(ii,1,ind1)==1,regnm(ii,2,ind1)==1) && or(regArcnm(ii,1,ind1)==1,regArcnm(ii,2,ind1)==1) ); %iii condition page 182; if there is one polarity region at either the 1st or 2nd thrshold AND there is one arc at either the 1st or 2nd thrshold
    cond2b1 = ( regnm(ii,1,ind2)==1 && regArcnm(ii,1,ind2)==1 ); %Condition ii part 1 from paper
    
    %cond2b2 = ( regnm(ii,1,ind2)==0 && regArcnm(ii,1,ind2)==0 && regY(ii,1,ind1)<=0 ); %Condition ii part 2 from paper
    cond2 = (cond2a && cond2b1);
    cond2b2 = ( regnm(ii,1,ind2)==0 && regArcnm(ii,1,ind2)==0 && regY(ii,1,ind1)<=0 ); %Condition ii part 2 from paper
    if flg==0
        cond2 = (cond2a && or(cond2b1,cond2b2));
    else
        cond2 = (cond2a && cond2b1);
    end
    cond3a = ( regcnt2org(ii,1)<0.3 ); %Condition 3: d <0.3 page 6 supp aEREMCOR
    cond3b = 0;
    if( cond3a==0 && cond2b==1 )
        reg12cnt2org(ii) = lne2pt([regX(ii,1,kk) regY(ii,1,kk) 0],[regX(ii,2,kkp) regY(ii,2,kkp) 0],[0 0 0]); %find distance between these points d'
        if( lt(reg12cnt2org(ii),0.3) ) % d'<0.3
            cond3b = 1;
        end
    end
    cond3 = ( cond3a || and(cond3b,condlgmopkspt(ii)) ); % save locs where either cond3a or (cond3b and the bcg identified ICs from spectrum analysis) are true
    cond4a = ( all(gt(regA(ii,1,1:2),0.1)) && all(gt(regArcA(ii,1,1:2),0.25)) ); %Condition 4: sets the minimum area of the polarity arc region parameters from table S4
    cond4b = ( max(regA(ii,1,1:2))>0.1 && min(regA(ii,1,1:2))>0.05 && max(regArcA(ii,1,1:2))>0.25 && min(regArcA(ii,1,1:2))>0.08 ); %Condition 4: sets the minimum area of the secondary polarity region parameters from table S4
    cond4c = 0;
    if( cond4a==0 && cond4b==0 && cond2b==1 ) % if all of condition 4 is met...
        if( regA(ii,1,kk)>0.1 && regA(ii,2,kkp)>0.1 && regArcA(ii,1,kk)>0.25 && regArcA(ii,2,kkp)>0.12 ) %and if the area of the regions are above a certain percentage
            cond4c = 1; % yes to Condition 4
        end
    end
    cond4 = ( regA(ii,2,ind1)>0.25 && regArcA(ii,2,ind1)>.025 ); %Condition 4 from Cardiac Paper is where these areas are above a certain percentage
    if( cond1 && cond2 && cond3 && cond4 ) %if all conditions are met
        if((cond3a||cond3b) && (cond4a||cond4b||cond4c)) %if one of the 3rd and 4th conditions are met
            condrmmot(ii) = 1; % note that location with a 1
        end
    end

    if flg==1
        if( (cond1 && cond2 && cond4) && ~cond5) %||cond7%if condition 1, condition 2, and condition 4 are met
            cbtpind(ii) = 1; % cbtpind tells us the IC numbers that are classified as BCG
        end
        if cond6
            Bipflg(ii)=1;
        end
    else
        if (cond1 && cond2 && cond4) %||cond7%if condition 1, condition 2, and condition 4 are met
            cbtpind(ii) = 1; % cbtpind tells us the IC numbers that are classified as BCG
        end
    end
end

% removal of blink and saccade components
condbk = zeros(1,icnum);
condsc = zeros(1,icnum);
condbk(intersect(find(condbkspt==1),find(condbktp==1))) = 1;
mofrqind(condbk==1) = 0;
motpind(condbk==1) = 0;
condlgmopkspt(condbk==1) = 0;
cbtpind(condbk==1) = 0;
condsc(intersect(find(condscspt==1),find(condsctp==1))) = 1;
mofrqind(condsc==1) = 0;
motpind(condsc==1) = 0;
condlgmopkspt(condsc==1) = 0;
cbtpind(condsc==1) = 0;

% signal contribution analysis
%rmpos
rmsignifind = zeros(1,icnum);
cbsignifind = zeros(1,icnum);
signifcond = zeros(icnum,nbch);
spkdurintv = 0.04/mriperd*szsm; %number of data point in 40ms duration
rmpos = zeros(szic(1),szic(2));
rmneg = zeros(szic(1),szic(2));
cbpos = zeros(szic(1),szic(2));
cbneg = zeros(szic(1),szic(2));
for ii=1:icnum
    rmpos(ii,ic(ii,1:szsm)>icave(ii)+4*icstd(ii)) = 1;
    rmneg(ii,ic(ii,1:szsm)<icave(ii)-4*icstd(ii)) = 1;
    cbpos(ii,ic(ii,1:szsm)>icave(ii)+0.1*icstd(ii)) = 1;
    cbneg(ii,ic(ii,1:szsm)<icave(ii)-0.1*icstd(ii)) = 1;
end
eachperd = floor(10*srte);
perdnm = floor(mriperd*srte/eachperd);
maxsig = zeros(nbch,1);
meansig = zeros(nbch,1);
stdsig = zeros(nbch,1);
minstdsig = zeros(nbch,1);
rmsig = zeros(nbch,2);
nrmperd = zeros(nbch,szsm);
for jj=1:nbch
    maxsig(jj) = max(abs(mixsig(jj,1:szsm)));
    meansig(jj) = mean(mixsig(jj,1:szsm));
    stdsig(jj) = std(mixsig(jj,1:szsm));
    minstdsig(jj) = std(mixsig(jj,1:eachperd));
    for kk=2:perdnm
        minstdsig(jj) = min(minstdsig(jj),std(mixsig(jj,(kk-1)*eachperd+1:kk*eachperd)));
    end
    rmsig(jj,1) = meansig(jj) + 4 * stdsig(jj);
    rmsig(jj,2) = meansig(jj) - 4* stdsig(jj);
    nrmperd(jj,mixsig(jj,1:szsm)<=meansig(jj)+4*minstdsig(jj) & mixsig(jj,1:szsm)>=meansig(jj)-4*minstdsig(jj) ) = 1;
end
flgSingleChArt=zeros(size(ic,1),1);
for ii = 1:icnum
    [possectnm,possectrge] = cont1(rmpos(ii,1:szsm));
    [negsectnm,negsectrge] = cont1(rmneg(ii,1:szsm));
    selic = zeros(size(ic));
    selic(ii,1:szsm) = ic(ii,1:szsm);
    corrsig = mixsig - A*selic;
    if flg==1
        subsig=A*selic;
        [spt1,frq1] = pwelchcal(subsig,srte,512,256); %spt is the power spectral density curve, frq is the freq. of the PSD curve
        spt1_sum = sum(spt1,2);
        [spt1sum spt1sumind] = sort(spt1_sum,'descend');
        spt1maxes = spt1sum(1:3);
        if (spt1maxes(1)>5*spt1maxes(2)||spt1maxes(1)>10*spt1maxes(3)) && ickurt(ii)>4  && alpharatio<max(max(deltaratio,thetaratio),betaratio)
            
            flgSingleChArt(ii)=1;
        end
    end
    %{
         spt1maxesind = spt1sumind(1:3);
         max3vals = spt1_sum(spt1maxesind(1:3)); %values from the top 3 max values
             for w = 1:size(spt1maxes,1)
                 %spt1ratios =
             end

    %}
    RapidMotionCondition = ( ickurt(ii)>5 && (sum(rmpos(ii,1:szsm))>0 || sum(rmneg(ii,1:szsm))>0) );
    if( RapidMotionCondition )
        cond1 = zeros(nbch,possectnm);
        for kk=1:possectnm
            if(possectrge(kk,2)-possectrge(kk,1)>spkdurintv) % check if the duration is more than 40 ms
                Diff = max(abs(mixsig(:,possectrge(kk,1):possectrge(kk,2))-corrsig(:,possectrge(kk,1):possectrge(kk,2))),[],2); %delta_jk; The signal reduction at the iith IC
                RelDiff = Diff ./ max(abs(mixsig(:,possectrge(kk,1):possectrge(kk,2))),[],2);%delta_jk/Vjk,max
                AbsDiff = Diff ./ max(maxsig(:),100*ones(size(maxsig))); %????
                Ave1 = abs(mean(mixsig(:,possectrge(kk,1):possectrge(kk,2)),2));
                Ave2 = abs(mean(corrsig(:,possectrge(kk,1):possectrge(kk,2)),2));
                SignalSpike = or(gt(max(mixsig(:,possectrge(kk,1):possectrge(kk,2)),[],2),rmsig(:,1)),lt(min(mixsig(:,possectrge(kk,1):possectrge(kk,2)),[],2),rmsig(:,2)));
                Significance = ( ( lt(Ave2,0.9*Ave1) & gt(RelDiff,diffr(1)*ones(size(RelDiff))) ) & ( gt(AbsDiff,diffr(2)*ones(size(AbsDiff))) | gt(Diff,200*ones(size(Diff))) ) );
                cond1(:,kk) = and( Significance(:), SignalSpike(:) );
                for jj=1:nbch
                    if(cond1(jj,kk)==1)
                    end
                end
            end
        end
        cond1(cond1(:, sum(cond1(:,:),1) < 2) == 1) = 0;
        signifcond(ii,:) = ( sum(cond1,2)>0 );
        
        cond1 = zeros(nbch,negsectnm);
        for kk=1:negsectnm
            if(negsectrge(kk,2)-negsectrge(kk,1)>spkdurintv)
                Diff = max(abs(mixsig(:,negsectrge(kk,1):negsectrge(kk,2))-corrsig(:,negsectrge(kk,1):negsectrge(kk,2))),[],2);
                RelDiff = Diff ./ max(abs(mixsig(:,negsectrge(kk,1):negsectrge(kk,2))),[],2);
                AbsDiff = Diff ./ max(maxsig(:),100*ones(size(maxsig)));
                Ave1 = abs(mean(mixsig(:,negsectrge(kk,1):negsectrge(kk,2)),2));
                Ave2 = abs(mean(corrsig(:,negsectrge(kk,1):negsectrge(kk,2)),2));
                SignalSpike = or(gt(max(mixsig(:,negsectrge(kk,1):negsectrge(kk,2)),[],2),rmsig(:,1)),lt(min(mixsig(:,negsectrge(kk,1):negsectrge(kk,2)),[],2),rmsig(:,2)));
                Significance = ( ( lt(Ave2,0.9*Ave1) & gt(RelDiff,diffr(1)*ones(size(RelDiff))) ) & ( gt(AbsDiff,diffr(2)*ones(size(AbsDiff))) | gt(Diff,200*ones(size(Diff))) ) );
                cond1(:,kk) = and( Significance(:), SignalSpike(:) );
                for jj=1:nbch
                    if(cond1(jj,kk)==1)
                    end
                end
            end
        end
        cond1(cond1(:, sum(cond1(:,:),1) < 2) == 1) = 0;
        signifcond(ii,:) = or( signifcond(ii,:), transpose(sum(cond1,2)>0) );
    end
    indrm = ( signifcond(ii,ignsignif==0)==1 );
    if( condrmmot(ii) && sum(indrm)>=6 )
        rmsignifind(ii) = 1;
    end
end

% bcg ic
cond3 = zeros(icnum,nbch);
for ii=1:icnum
    selic = zeros(size(ic));
    selic(ii,1:szsm) = ic(ii,1:szsm);
    corrsig = mixsig - A*selic; % Corrected Signal after removing IC ii
    if( sum(cbpos(ii,1:szsm))>0 && sum(cbneg(ii,1:szsm))>0 )
        for jj=1:nbch
            posave1 = mean(abs( mixsig(jj,cbpos(ii,1:szsm)==1 & nrmperd(jj,1:szsm)==1) ));
            posave2 = mean(abs( corrsig(jj,cbpos(ii,1:szsm)==1 & nrmperd(jj,1:szsm)==1) ));
            negave1 = mean(abs( mixsig(jj,cbneg(ii,1:szsm)==1 & nrmperd(jj,1:szsm)==1) ));
            negave2 = mean(abs( corrsig(jj,cbneg(ii,1:szsm)==1 & nrmperd(jj,1:szsm)==1) ));
            cond2 = ( min(posave2/posave1,negave2/negave1)<0.95 && 0.5*(posave2/posave1+negave2/negave1)<0.97 );
            if(cond2)
                signifcond(ii,jj) = signifcond(ii,jj) + 2;
            end
            cond3(ii,jj) = ( 0.5*(posave2/posave1+negave2/negave1)<0.8 );
        end
        indcb = ( signifcond(ii,:)==2 | signifcond(ii,:)==3 );
        if( condlgmopkspt(ii) && sum(indcb)>=1 )
            cbsignifind(ii) = 1;
        end
    end
    signifcond(ii,find(cond3(ii,:)==1)) = signifcond(ii,find(cond3(ii,:)==1)) + 4;
    indrsp = ( signifcond(ii,ignsignif==0)>=4 );
end

% selection combination
cbicind = ( condlgmopkspt==1 & cbtpind==1 & cbsignifind==1 )|(Bipflg==1 & alphababy == 0);
if flg==0
    saccade=[];
    blink=[];
    topomap=[];
    spectrumy=[];
    tpblink=[];
    tpsac=[];
    smolregion=[];
    singchan=[];
    muscleloc=[];
else
    topomap = cbtpind;
    spectrumy = condlgmopkspt;
    blink = condbk;
    saccade = condsc;
    smolregion = toosmall;
    tpblink = condbktp;
    tpsac = condsctp;
    singchan=flgSingleChArt;
    muscleloc=gammaflg;
end
end
