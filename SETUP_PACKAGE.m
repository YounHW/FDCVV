function [root,upstream]=SETUP_PACKAGE(method,upstream)
% Set the original framework location once per fresh MATLAB session.
root=fileparts(mfilename('fullpath'));
if nargin<2||isempty(upstream),upstream=getenv('FDCVV_UPSTREAM_ROOT');end
if isempty(upstream),upstream=fullfile(fileparts(root),'StateOfTheArtAdaptiveSampling');end
assert(isfolder(fullfile(upstream,'src')),'Set FDCVV_UPSTREAM_ROOT to the pinned original repository folder. See README.');
upstream=char(java.io.File(upstream).getCanonicalPath());
lock=jsondecode(fileread(fullfile(root,'upstream_lock.json')));
for i=1:numel(lock.files)
 p=fullfile(upstream,lock.files(i).file);assert(isfile(p),'Missing dependency: %s',p);
 txt=fileread(p);txt=strrep(txt,sprintf('\r\n'),sprintf('\n'));
 if ~isempty(txt)&&txt(1)==char(65279),txt=txt(2:end);end
 md=java.security.MessageDigest.getInstance('SHA-256');md.update(unicode2native(txt,'UTF-8'));
 h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
 assert(strcmp(h,lock.files(i).sha256_lf),'Dependency differs from pinned revision: %s',p);
end
setenv('FDCVV_UPSTREAM_ROOT',upstream);
restoredefaultpath;pp=strsplit(path,pathsep);
pp=pp(cellfun(@(p)startsWith(lower(p),lower(matlabroot)),pp));path(strjoin(pp,pathsep));
addpath(root);addpath(genpath(fullfile(upstream,'src')),'-begin');
addpath(genpath(fullfile(root,'matlab','base')),'-begin');
if ismember(method,{'FDCVV','LC2015'}),addpath(fullfile(root,'matlab','lc2015'),'-begin');end
if strcmp(method,'FDCVV'),addpath(fullfile(root,'matlab','final_fdcvv'),'-begin');end
global plotFlag isFEAFlag;plotFlag=0;isFEAFlag=false;maxNumCompThreads(1);
expected=fullfile(root,'matlab','base','OK_model.m');
if ismember(method,{'FDCVV','LC2015'}),expected=fullfile(root,'matlab','lc2015','OK_model.m');end
assert(strcmp(which('OK_model'),expected),'Use a fresh MATLAB session');
end
