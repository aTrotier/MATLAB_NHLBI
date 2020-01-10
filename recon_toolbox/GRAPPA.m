function [recon, kspace_filled] = GRAPPA(kspace, mask, af)
% [recon, sigrecon] = GRAPPA(kspace, mask, af)
% Wrapper for DNF's version of M Griswold's grappa function
% integrated implementation 

% figure, imagesc(mask);
mask = mask(1,:);

acs = kspace(:, mask==3 ,:);
kspace = kspace(:,1:af:end,:);
sig = permute(kspace,[3 2 1]);
acs = permute(acs,[3 2 1]);

[recon, kspace_filled]=GRAPPA_open(sig,acs,af);
recon = permute(recon,[3 2 1]);
kspace_filled = permute(kspace_filled,[3 2 1]);

% montage_RR(abs(recon ))
% montage_RR(abs(kspace_filled ))

end

function [recon, sigrecon]=GRAPPA_open(sig,acs,af);
% function [recon, sigrecon]=GRAPPA_open(sig,acs,af);
%
%   Please read the license text at the bottom of this program. By using this program, you
%   implicity agree with the license.
%
%   The main points of the license:
%   1) This code is strictly for non-commercial applications. The code is
%   protected by
%      multiple patents.
%   2) This code is strictly for research purposes, and should not be used in any
%      diagnostic setting.
%
%   22.10.2004  Mark Griswold (mark@physik.uni-wuerzburg.de)


[nc,ny,nx]=size(sig);
[nc,nyacs,nxacs]=size(acs);

src=zeros((nyacs-af)*(nxacs-2),nc*6);
targ=zeros((nyacs-af)*(nxacs-2),nc*(af-1));

cnt=0;
for xind=2:nxacs-1,
      for yind=1:nyacs-af,

        cnt=cnt+1;
        src(cnt,:)=reshape(acs(:,yind:af:yind+af,xind-1:xind+1),1,nc*6);
        targ(cnt,:)=reshape(acs(:,yind+1:yind+af-1,xind),1,nc*(af-1));

    end
end

ws=pinv(src)*targ;

sigrecon=zeros(nc,ny*af,nx);
sigrecon(:,1:af:end,:)=sig;

for xind=2:nx-1,
    for yind=1:af:ny*af-af,

        srcr=reshape(sigrecon(:,yind:af:yind+af,xind-1:xind+1),1,nc*6);
        sigrecon(:,yind+1:yind+af-1,xind)=reshape(srcr*squeeze(ws),[nc af-1]);

    end
end

recon=circshift(ifft(ifft(circshift(sigrecon,[0 round(af*ny./2) round(nx./2)]),[],2),[],3),[0 round(af*ny./2) round(nx./2)]);

end

%% license Agreement
% You should carefully read the following terms and conditions before installing or using the
% software. Unless you have entered into a separate written license agreement with
% Universit?t W?rzburg providing otherwise, installation or use of the software indicates your
% agreement to be bound by these terms and conditions.
%
% Use of the software provided with this agreement constitutes your acceptance of these terms.
% If you do NOT agree to the terms of this agreement, promptly remove the software together
% with all copies from your computer. User's use of this software is conditioned upon compliance
% by user with the terms of this agreement.
%
% Upon ordering, downloading, copying, installing or unencrypting any version of the software, you
% are reaffirming that you agree to be bound by the terms of this agreement.
%
% License to use
%
% Universit?t W?rzburg grants to you a limited, non-exclusive, non-transferable and non-assignable
% license to install and use this software for research purposes. Use of this software for any
% diagnostic imaging procedure is strictly forbidden.
%
% License to distribute
%
% Please feel free to offer the non-commercial version of this software on any website, CD, or
% bulletin board, demonstrate the non-commercial version of the software and its capabilities, or
% give copies of the non-commercial version of the software to other potential users, so that others
% may have the opportunity to obtain a copy for use in accordance with the license terms contained
% here.
%
% You agree you will only copy the non-commercial version of the software in whole with this
% license and all delivered files, but not in part.
%
% Termination
%
% This license is effective until terminated. You may terminate it at any point by destroying
% the software together with all copies of the software.
%
% If you have acquired a non-commercial version, the license granted herein shall automatically
% terminate if you fail to comply with any term or condition of this Agreement.
%
% Also, Universit?t W?rzburg has the option to terminate any license granted herein if you fail
% to comply with any term or condition of this Agreement.
%
% You agree upon such termination to destroy the software together with all copies of the software.
%
%
% Copyright
%
% The software is protected by copyright law. You acknowledge that no title to the intellectual
% property in the software is transferred to you. You further acknowledge that title and full
% ownership rights to the software will remain the exclusive property of Universit?t W?rzburg,
% and you will not acquire any rights to the software except as expressly set forth in this
% license. You agree that any copies of the software will contain the same proprietary notices
% which appear on and in the software.
%
% Rent, lease, loan
%
% You may NOT rent, lease or loan the software without first negotiating a specific license
% for that purpose with Universit?t W?rzburg.
%
% No warranties
%
% Universit?t W?rzburg does NOT warrant that the software is error free. Universit?t W?rzburg
% disclaims all warranties with respect to the software, either express or implied, including
% but not limited to implied warranties of merchantability, fitness for a particular purpose and
% noninfringement of third party rights. The software is provided "AS IS."
%
% No liability for consequential damages
%
% In no event will Universit?t W?rzburg be liable for any loss of profits, business, use, or data
% or for any consequential, special, incidental or indirect damages of any kind arising out of
% the delivery or performance or as a result of using or modifying the software, even if
% Universit?t W?rzburg has been advised of the possibility of such damages. In no event will
% Universit?t W?rzburg's liability for any claim, whether in contract, negligence, tort or any
% other theory of liability, exceed the license fee paid by you, if any.
% The licensed software is not designed for use in high-risk activities requiring fail-safe
% performance. Universit?t W?rzburg disclaims any express or implied warranty of fitness for
% high-risk activities.
%
% Severability
%
% In the event of invalidity of any provision of this license, the parties agree that such
% invalidity shall not affect the validity of the remaining portions of this license.
%
