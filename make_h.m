% make_h_simple
% Takes a scalar set of parameters and returns an array with values for h.
% h(X,a,i,j,k) is the susceptibility of individuals to strain X dependent on their
% age and their antibody levels. Values are referenced to h(X,1,1,1,1)=1, i.e.
% those in the youngest age group with no detectable titres always have a relative
% suceptibility of 1
% Written by Hsiang-Yu Yuan (sean.yuan@cityu.edu.hk)
function [ h ] = make_h( pa)

maxa = pa.maxa; %number of age groups
maxi = pa.maxi; %Ab level to strain A
maxj = pa.maxj; %Ab level to strain B
maxk = pa.maxk; %Ab level to strain C
ka = pa.ka; %shape parameter
kb = pa.kb; %0.5 protection
h = zeros(maxa,maxi,maxj,maxk);

%option 1 from 100% to 0%
%before = 1./(1+exp(ka*[[1:maxi]-kb]))-1./(1+exp(ka*[maxi-kb]));
%ratio = 1/(1./(1+exp(ka*[1-kb]))-1./(1+exp(ka*[maxi-kb])));
%h = before.*ratio; 

%option 2 from 100% to lower bound
before = 1./(1+exp(ka*[[1:maxi]-kb]))-1./(1+exp(ka*[maxi-kb]));
ratio = (1-1./(1+exp(ka*[maxi-kb])))/(1./(1+exp(ka*[1-kb]))-1./(1+exp(ka*[maxi-kb])));
h = before.*ratio + 1./(1+exp(ka*[maxi-kb])); 

%option 3
%    for a=1:maxa
%        for i=1:maxi
%            for j=1:maxj
%                for k=1:maxk
%                        % h is linearly decreasing
%                        %Change from 1.3 to 1.25
%                        h(a,i,j,k)  = 1./(1+exp(ka*[i-kb])); %titre 1:25.6 represent titre level 3.3561
%                end
%            end
%        end
%    end
end
