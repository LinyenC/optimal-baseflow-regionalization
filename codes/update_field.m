% Yiwen Mei (yiwen.mei@uconn.edu)
% CIRCA, University of Connecticut
% First version on 5/1/2022
% Last updated on 10/18/2022

%% Functionality:
% This code fills the parameter values marked as NaN for the C_md or BF_md table
%  with the calibrated values.

%% Input
% F_md: the C_md or the BF_md table;
%        - C_md for CFit has three parameters, "c_str", "c_f", and "n";
%        - BF_md for RDF has one parameter, "BFI_x";
%        - BF_md for SMM has two parameters, "N_d" and "rt";
%        - BF_md for FSMM has three parameters, "BFI_x", "N_d", and "rt";
%        - BF_md for BRM has two parameters, "f_bmp", and "k_rise".
% Prms: calibrated parameters used to replace the NaN in F_md (it can be supplied
%         as a 1-by-n table with the n calibrated parameters matching the number
%         of NaN field in F_md; it can be supplied as a structure with the same
%         field as F_md without NaN);

%% Output
% F_md_o: updated F_md with the NaN field filled;
%  Fdn  : the names of the field got updated;
%  Val  : the values of the field got updated.

function [F_md_o,Fdn,Val]=update_field(F_md,Prm)
%% Check the inputs
narginchk(2,2);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'F_md',@(x) validateattributes(x,{'table'},{'nonempty'},mfilename,'F_md'));
if istable(Prm)
  addRequired(ips,'Prm',@(x) validateattributes(x,{'table'},{'nonempty'},mfilename,'Prm'));
  fld_nm=F_md.Properties.VariableNames;
  fld_nm_o=Prm.Properties.VariableNames;
  if ~isequal(fld_nm_o,fld_nm)
    error('Field names of F_md and Prm must be identical');
  end
elseif isvector(Prm)
  addRequired(ips,'Prm',@(x) validateattributes(x,{'double'},{'nonempty','nonnan'},mfilename,'Prm'));
  k=table2array(varfun(@isnumeric,F_md));
  nV_cal=sum(isnan(table2array(F_md(:,k))));
  if nV_cal~=length(Prm)
    error('Number of variable to update in F_md does not match the number of value supplied by Prm');
  end
end

parse(ips,F_md,Prm);
clear ips

%% Fill the calibrated parameters
F_md_o=F_md;
fld_nm=F_md.Properties.VariableNames;
Val=[];
Fdn=[];
a=1;
for j=1:length(fld_nm)
  f_nm=F_md{1,j};

  if isnumeric(f_nm)
    if isnan(f_nm)
      if istable(Prm)
        val=Prm{1,j};
      else
        val=Prm(a);
      end
      Val=[Val val];
      Fdn=[Fdn fld_nm{j}];
      F_md_o{1,j}=val;
      a=a+1;
    end
  end
end
end
