function [oriabu,bioabu,oriiso,bioiso,biopart1] = iturbo2_PETM690(abu,iso,mxl,numb)
% The MATLAB program iTURBO2 can be used to simulate the effects of
% bioturbation on single sediments particles.
%
% HERE adapted for PETM experiments to pass back biopart1
%
% ABU = series of abundances of carrier type 1 down core
% MXL = series of mixed layer thicknesses down core
% NUMB = number of carriers to be measured
% ORIABU = original abundances of both carrier types 1 and 2
% BIOABU = bioturbated abundances of both carriers types 1 and 2
% ORIISO = original isotope signature of both carrier  types 1 and 2
% BIOISO = bioturbated isotope signature of both carrier types 1 and 2
% BIOPART1 = values of single foram particles that get measured

% mainly from Martin Trauth 18 July 2012
nrows = max(mxl); % +10;        % beginning of artifical sediment column (will grow each time-step)
ncols = max(abu); % +50;        % max. # of foram tests in abu + 50

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
%
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
    bioiso(i,1) = nanmean(biopart1(i,:));   % calculate bioturbated isotope mean of the numb (50) values
    bioiso(i,2) = nanmean(biopart2(i,:));   % calculate bioturbated isotope mean of the numb (50) values
end
% flip abu and iso upside down, otherwise sediment column grows downward
oriabu = flipud(oriabu);
oriiso = flipud(oriiso);
bioabu = flipud(bioabu);
bioiso = flipud(bioiso);
biopart1 = flipud(biopart1);







