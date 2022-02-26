function [oriabu,bioabu,oriiso,bioiso] = iturbo2_plus_TM(abu,iso,mxl,numb, TM)
% The MATLAB program iTURBO2 (re-) introduced the use of a transtion matrix (TM)
% and can be used to simulate the effects of bioturbation on single sediments particles.
%
% ABU = series of abundances of carrier type 1 down core
% MXL = series of mixed layer thicknesses down core
% NUMB = number of carriers to be measured
% TM = logical TRUE: use transition matrix; FALSE: use turbo2 homogenous mixing
% ORIABU = original abundances of both carrier types 1 and 2
% BIOABU = bioturbated abundances of both carriers types 1 and 2
% ORIISO = original isotope signature of both carrier  types 1 and 2
% BIOISO = bioturbated isotope signature of both carrier types 1 and 2

if(~TM)
    % From Martin Trauth 18 July 2012
    nrows = max(mxl);%+10;        % beginning of artifical sediment column (will grow each time-step)
    ncols = max(abu)+50;        % max. # of foram tests in abu + 50
    
    sedabu = NaN(nrows,ncols);
    sediso = NaN(nrows,ncols);
    % h = waitbar(0,'Mixing Process ...');
    for i = 1:length(abu)           % mixing for loop
        %  waitbar(i/length(abu))
        rncols = randperm(ncols);     % random permutation of 1, ..., ncols
        sedabu(size(sedabu,1)+1,1:ncols) = 2*ones(1,ncols);  % new layer of on top with 2.0
        sedabu(size(sedabu,1),1:abu(i)) = ones(1,abu(i));     % make 1.0 for abu
        sedabu(size(sedabu,1),:) = sedabu(size(sedabu,1),rncols); % mix the new layer horizontally
        sediso(size(sediso,1)+1,1:ncols) = iso(i)*ones(1,ncols); %  new layer of isotope values (both species with same value)
        for j = 1 : ncols     % actual mixing process
            z = randperm(mxl(i));     % random permutation of 1, ..., xml
            ns(1:mxl(i),j) = sedabu(size(sedabu,1)-mxl(i)+1:size(sedabu,1),j);        % copy mixed layer of column j
            sedabu(size(sedabu,1)-mxl(i)+1:size(sedabu,1),j) = ns(z,j);       % mix copy vertically within mixed layer
            ni(1:mxl(i),j) = sediso(size(sediso,1)-mxl(i)+1:size(sediso,1),j);    % same mixing with isotopes
            sediso(size(sediso,1)-mxl(i)+1:size(sediso,1),j) = ni(z,j);       % same mixing with isotopes
        end
    end
    % close(h)
    
else % USE TRANSITION MATRIX
    %%
    
    %    sedrate = 1e-4; % cm/yr
    %    dt = 1./sedrate; % yr
    dt = max(abu);   % was: abu(1);
    dtint = round(dt);
    dx = 1.; % cm
    % read transition matrix
    mxltxt = num2str(mxl(1));
    mxl_filename = ['transmtx/trans-equal_local_',mxltxt ,'cm.txt'];
    fileID=fopen(mxl_filename,'r');
%     % was before:
%     if(mxl(1)==5)
%         fileID=fopen('transmtx/trans-equal_local_5cm.txt','r');
%     elseif(mxl(1)==10)
%         fileID=fopen('transmtx/trans-equal_local_10cm.txt','r');
%     elseif(mxl(1)==20)
%         fileID=fopen('transmtx/trans-equal_local_20cm.txt','r');
%     end
    sizeA=[mxl(1) mxl(1)];
    formatSpec = '%f';
    A = fscanf(fileID,formatSpec,sizeA);
    A = A';
    A = A*dt;
    fclose(fileID);
    % end reading
    exrows = 0;
    excols = 0;
    % mxl(:) = 5; % fixing mxl
    nrows = max(mxl)+exrows;
    ncols = max(abu)+excols;
    % ncols = 1500; % fixing ncols
    sedabu = NaN(nrows,ncols);
    sediso = NaN(nrows,ncols);
    %h = waitbar(0,'Mixing Process ...');
    for i = 1:length(abu) %  length(abu) is number of rows in abu;
        %  the same as the other input variables;
        %  this is the number of layers deposited
        %  waitbar(i/length(abu))
        sedabu(size(sedabu,1)+1,1:ncols) = 2*ones(1,ncols);
        % size(sedabu,1)+1; this exceed the size of sedabu;
        % but this add 2 to the (added)bottom of sedabu
        % meaning all species 2 at the bottom of sedabu?
        sedabu(size(sedabu,1),1:abu(i)) = ones(1,abu(i));
        % fill 1st to i-th column with abu(i) individuals of 1 (species 1) at the bottom
        % the rest remains 2 (species 2)
        sediso(size(sediso,1)+1,1:ncols) = iso(i)*ones(1,ncols);
        % adding one layer to the bottom of sediso; assigned all isotope value of
        % species 1 during i-th deposition
        
        ns(1:mxl(i),:) = sedabu(size(sedabu,1)-mxl(i)+1:size(sedabu,1),:); % creating a mirror of mix layer
        ni(1:mxl(i),:) = sediso(size(sediso,1)-mxl(i)+1:size(sediso,1),:);
        for k = 1:mxl(i)
            for j=1:mxl(i)
                if j==k
                    continue
                end
                if round(A(j,k))==0
                    continue
                end
                
                z = randperm(ncols,round(A(j,k))); % round(A(j,k)) represents the number of partilces exchanged
                % between layers j and k
                % columns of mixing is
                % randomly chosen out of ncols
                sedabu(size(sedabu,1)+1-j,z) = ns(size(ns,1)+1-k,z);
                sedabu(size(sedabu,1)+1-k,z) = ns(size(ns,1)+1-j,z);
                % random mixing in vertical direction
                sediso(size(sediso,1)+1-j,z) = ni(size(ni,1)+1-k,z);
                sediso(size(sediso,1)+1-k,z) = ni(size(ni,1)+1-j,z);
                % the same is done for sediso for vertical mixing
                clear z
            end
        end
        clear ns ni
        
        %   end
    end
    
end
% flipud(sedabu(:,1:5))
% flipud(sediso(:,1:5))
%
clear ns ni i j mxl z
sedabu = sedabu(nrows+1:end,:); % remove initial sediment pile
sediso = sediso(nrows+1:end,:); % remove initial sediment pile

% calculate abundances of Species 1 & 2 in each layer
oriabu(:,1) = abu;          % orignal abundance of Sp. 1 before mixing
oriabu(:,2) = ncols-abu;    % orignal abundance of Sp. 2 before mixing
bioabu(:,1) = sum(sedabu==1,2); % abundance of Sp. 1 after mixing/bioturbation
bioabu(:,2) = sum(sedabu==2,2); % abundance of Sp. 2 after mixing/bioturbation

% calculate isotope values of Sp 1 and 2
oriiso(:,1) = iso;  % original signature
oriiso(:,2) = iso;  % original signature
bioiso1 = sediso;   % copy of mixed isotopes
bioiso2 = sediso;   % copy of mixed isotopes
bioiso1(sedabu~=1) = NaN;   % delete particles which are not 1
bioiso2(sedabu~=2) = NaN;   % delete particles which are not 2

biopart1 = NaN(size(bioiso1));
biopart2 = NaN(size(bioiso2));
for i = 1:length(abu)
    biopart1(i,1:bioabu(i,1)) = bioiso1(i,isnan(bioiso1(i,:))==0);  % replace by non-NAN elements and have values in first abu elements
    biopart2(i,1:bioabu(i,2)) = bioiso2(i,isnan(bioiso2(i,:))==0);  % replace by non-NAN elements and have values in first abu elements
end
biopart1 = biopart1(:,1:numb);  % reduce array to number of carriers to be picked/measured
biopart2 = biopart2(:,1:numb);  % reduce array to number of carriers to be picked/measured
for i = 1:length(abu)
    bioiso(i,1) = mean(biopart1(i,:), 'omitnan');   % calculate bioturbated isotope mean of the numb (50) values
    bioiso(i,2) = mean(biopart2(i,:), 'omitnan');   % calculate bioturbated isotope mean of the numb (50) values
end
% flip abu and iso upside down, otherwise sediment column grows downward
oriabu = flipud(oriabu);
oriiso = flipud(oriiso);
bioabu = flipud(bioabu);
bioiso = flipud(bioiso);

