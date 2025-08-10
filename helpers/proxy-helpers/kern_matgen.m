function A = kern_matgen(i,j,srcinfo,targinfo,coefs,kernfun,idspmat)
%
% 
%  kern_matgen
%    This subroutine returns the matrix block A(I,J) corresponding 
%    to the discretization points indexed by i\in I and j\in J. 
%  
%  Syntax: 
%    A = kern_matgen(i,j,srcinfo,targinfo,kern)
%
%  Input arguments:
%    * i: row indices
%    * j: column indices
%    * srcinfo: sources
%    * targino: tarks
%    * kern: kern(s,t) evaluates at sources and targets
%     
     if isempty(i) || isempty(j)
        A = zeros(length(i),length(j));
        return;
     end
     srcuse = [];
     srcuse.r = srcinfo.r(:,j);
     if isfield(srcinfo,'n')
        srcuse.n = srcinfo.n(:,j);
     end
     if isfield(srcinfo,'d')
        srcuse.d = srcinfo.d(:,j);
     end
     if isfield(srcinfo,'d2')
        srcuse.d2 = srcinfo.d2(:,j);
     end
      if isfield(srcinfo,'wts')
        srcuse.wts = srcinfo.wts(j);
     end
     targuse = [];
     targuse.r = targinfo.r(:,i);

     A1 = kernfun(srcuse, targuse);
     A = sum(coefs(i,:,:) .* A1,3);
     A = bsxfun(@times,A,srcinfo.wts(j).');
     A = A + idspmat(i,j);
     
end
