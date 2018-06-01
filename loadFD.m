%% available finite differences schemes
% L. LAURENT --  13/04/2018 -- luc.laurent@lecnam.net

% sources could be found here:
% https://bitbucket.com/luclaurent/gradFD
% https://github.com/luclaurent/gradFD
%

% gradFD - A toolbox to compute derivatives and hessians using finite differences
% Copyright (C) 2018  Luc LAURENT <luc.laurent@lecnam.net>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function R=loadFD(typeIn)
%from http://web.media.mit.edu/~crtaylor/calculator.html
%available technics and parameters
listFD.FD1.grad.steps=[0 1];
listFD.FD1.grad.coef=[-1 1];
listFD.FD1.grad.div=1;
listFD.FD1.hess.steps=[0 1 2];
listFD.FD1.hess.coef=[1 -2 1];
listFD.FD1.hess.div=1;
%
listFD.BD1.grad.steps=[-1 0];
listFD.BD1.grad.coef=[-1 1];
listFD.BD1.grad.div=1;
listFD.BD1.hess.steps=[-2 -1 0];
listFD.BD1.hess.coef=[1 -2 1];
listFD.BD1.hess.div=1;
%
listFD.FD2.grad.steps=[0 1 2];
listFD.FD2.grad.coef=[-3 4 -1];
listFD.FD2.grad.div=2;
listFD.FD2.hess.steps=[0 1 2 3];
listFD.FD2.hess.coef=[2 -5 4 -1];
listFD.FD2.hess.div=1;
%
listFD.BD2.grad.steps=[-2 -1 0];
listFD.BD2.grad.coef=[1 -4 3];
listFD.BD2.grad.div=2;
listFD.BD2.hess.steps=[-3 -2 -1 0];
listFD.BD2.hess.coef=[-1 4 -5 2];
listFD.BD2.hess.div=1;
%
listFD.CD2.grad.steps=[-1 1];
listFD.CD2.grad.coef=[-1 1];
listFD.CD2.grad.div=2;
listFD.CD2.hess.steps=[-1 0 1];
listFD.CD2.hess.coef=[1 -2 1];
listFD.CD2.hess.div=1;
%
listFD.FD3.grad.steps=[0 1 2 3];
listFD.FD3.grad.coef=[-11 18 -9 2];
listFD.FD3.grad.div=6;
listFD.FD3.hess.steps=[0 1 2 3 4];
listFD.FD3.hess.coef=[35 -104 114 -56 11];
listFD.FD3.hess.div=12;
%
listFD.BD3.grad.steps=[-3 -2 -1 0];
listFD.BD3.grad.coef=[-2 9 -18 11];
listFD.BD3.grad.div=6;
listFD.BD3.hess.steps=[-4 -3 -2 -1 0];
listFD.BD3.hess.coef=[-35 104 -114 56 -11];
listFD.BD3.hess.div=12;
%
listFD.FD4m.grad.steps=[-1 0 1 2 3];
listFD.FD4m.grad.coef=[-3 -10 18 -6 1];
listFD.FD4m.grad.div=12;
listFD.FD4m.hess.steps=[-1 0 1 2 3];
listFD.FD4m.hess.coef=[11 -20 6 4 -1];
listFD.FD4m.hess.div=12;
%
listFD.FD4.grad.steps=[0 1 2 3 4];
listFD.FD4.grad.coef=[-25 48 -36 16 -3];
listFD.FD4.grad.div=12;
listFD.FD4.hess.steps=[0 1 2 3 4 5];
listFD.FD4.hess.coef=[NaN];
listFD.FD4.hess.div=NaN;
%
listFD.BD4p.grad.steps=[-3 -2 -1 0 1];
listFD.BD4p.grad.coef=[-1 6 -18 10 3];
listFD.BD4p.grad.div=12;
listFD.BD4p.hess.steps=[-3 -2 -1 0 1];
listFD.BD4p.hess.coef=[-1 4 6 -20 11];
listFD.BD4p.hess.div=12;
%
listFD.BD4.grad.steps=[-4 -3 -2 -1 0];
listFD.BD4.grad.coef=[3 -16 36 -48 25];
listFD.BD4.grad.div=12;
listFD.BD4.hess.steps=[-5 -4 -3 -2 -1 0 1];
listFD.BD4.hess.coef=[NaN];
listFD.BD4.hess.div=NaN;
%
listFD.CD4.grad.steps=[-2 -1 1 2];
listFD.CD4.grad.coef=[1 -8 8 -1];
listFD.CD4.grad.div=12;
listFD.CD4.hess.steps=[-2 -1 0 1 2];
listFD.CD4.hess.coef=[-1 16 -30 16 -1];
listFD.CD4.hess.div=12;
%
listFD.FD5.grad.steps=[0 1 2 3 4 5];
listFD.FD5.grad.coef=[-137 300 -300 200 -75 12];
listFD.FD5.grad.div=60;
listFD.FD5.hess.steps=[0 1 2 3 4 5 6];
listFD.FD5.hess.coef=[NaN];
listFD.FD5.hess.div=NaN;
%
listFD.BD5.grad.steps=[-5 -4 -3 -2 -1 0];
listFD.BD5.grad.coef=[-12 75 -200 300 -300 137];
listFD.BD5.grad.div=60;
listFD.BD5.hess.steps=[-5 -4 -3 -2 -1 0 1];
listFD.BD5.hess.coef=[NaN];
listFD.BD5.hess.div=NaN;
%
listFD.CD6.grad.steps=[-3 -2 -1 1 2 3];
listFD.CD6.grad.coef=[-1 9 -45 45 -9 1];
listFD.CD6.grad.div=60;
listFD.CD6.hess.steps=[-3 -2 -1 0 1 2 3];
listFD.CD6.hess.coef=[2 -27 270 -490 270 -27 2];%fix
listFD.CD6.hess.div=180;
%
listFD.CD8.grad.steps=[-4 -3 -2 -1 1 2 3 4];
listFD.CD8.grad.coef=[3 -32 168 -672 672 -168 32 -3];
listFD.CD8.grad.div=840;
listFD.CD8.hess.steps=[-2 -1 0 1 2]; %fix
listFD.CD8.hess.coef=[-1 16 -30 16 -1];%fix
listFD.CD8.hess.div=12;%fix
%
if nargin==0
    typeIn=[];
end
%
if ~isempty(typeIn)
    listT=fieldnames(listFD);
    if ismember(typeIn,listT)
        R=listFD.(typeIn);
    else
        fprintf('Wrong type of finite differences\n');
        fprintf('Available technics:\n');
        cellfun(@(X)fprintf('%s\n',X),listT);
    end
else
    R=listFD;
end
end