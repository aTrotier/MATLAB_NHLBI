function [cg_cresult, cg_resid] = cg_RR(data, st, I, D, csm, weights, iter)
% Just spits back the end Magnitude image
%%
if nargin < 7
    iter = 25;
end

%% Initial Guess

%  a1 = I(:).*ismrm_encoding_non_cartesian_SENSE(D(:).*data(:),csm,st,weights,'transp');

% explicit
scale = sqrt(prod(prod(st.Kd))/numel(weights(:)));
samples = size(st.om,1);
coils = numel(data)/samples; 
data = reshape(data,samples,coils);

outp = (nufft_adj(D.*data, st))./(sqrt(prod(st.Kd)))*scale;
a1 = I.*sum(conj(csm) .* outp,3);

% % % a1 = ismrm_encoding_non_cartesian_SENSE(data(:),csm,st,weights,'transp');

%% CG LOOP INITIALISE
alpha = a1(:);
eps = 1e-5; % required accuracy
clear bapprox; bapprox{1} = 0;
p = a1(:);
r = a1(:);

%% LOOP

for i = 1:iter
    delta = (r' * r)/(alpha' * alpha);
    cg_resid(i) = delta;
    
    if delta < eps
        break;
    else
% % % %         q  = I*E'*D*E*I*p:
%         temp1 = D(:).*ismrm_encoding_non_cartesian_SENSE(I(:).*p(:),csm,st,weights,'notransp');
%         q = I(:).*ismrm_encoding_non_cartesian_SENSE(temp1,csm,st,weights,'transp');

    outp = repmat(reshape(I(:).*p, size(csm,1),size(csm,2)),[1 1 size(csm,3)]) .* csm;
    outp = (nufft(outp,st)./(sqrt(prod(st.Kd))))*scale;
    outp = (nufft_adj(outp .* repmat(weights,[1 coils]),st)./(sqrt(prod(st.Kd))))*scale;
    outp = I.*sum(conj(csm) .* outp,3);
    q = outp(:);

 
        
% % %          temp1 = ismrm_encoding_non_cartesian_SENSE(p(:),csm,st,weights,'notransp');
% % %          q = ismrm_encoding_non_cartesian_SENSE(temp1,csm,st,weights,'transp');
         
        % figure,imshow(abs([reshape(p, 256,256)/max(p) reshape(q, 256,256)/max(q)]),[]);
        
        bapprox{i+1} = bapprox{i} + ((r' * r)/(p' * q))*p;
        
        % figure('Name',['iter ' num2str(i)]),imshow((reshape(bapprox{i+1}.*conj(bapprox{i+1}), 128,128)),[]);
        
        rnew = r - ((r' * r)/(p' * q)) * q;
        p = rnew + ((rnew' * rnew)/(r' * r)) * p;
        r = rnew;
    end
end

%%

cg_cresult = reshape(bapprox{end}, floor(sqrt(length(bapprox{end}))), floor(sqrt(length(bapprox{end}))));

end