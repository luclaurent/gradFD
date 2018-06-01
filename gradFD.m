%% gradFD class for computing gradients using finite differences
% L. LAURENT --  12/12/2017 -- luc.laurent@lecnam.net

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


classdef gradFD < handle
    
    properties
        type='FD1';         % type of finite difference
        Xref;               % sample points on which the gradients will be calculated
        XevalG;             % points used for computing FD
        ZevalG;             % responses used for computing FD
        GZeval;             % values of gradients
        XevalH;             % points used for computing hessians with FD
        ZevalH;             % responses used for computing hessians with FD
        HZeval;             % values of hessians
        fun;                % considered function
        stepsDiff=1e-4;     % steps of FD
        dim;                % dimension (number of design variables)
        nS;                 % number of sample points
        nX;                 % number of points used for compute FD (duplicate points are removed)
    end
    properties (Access=private)
        confFD;             % configuration of the chosen scheme used for FD
        dupX;               % duplicate coordinates of Xeval (for reducing the number of evaluations of the function
        stepsDiffInternal;  % steps of FD
        ZevalGext;          % array of responses used in the case of external function
    end
    
    methods
        %% constructor of finite differences class
        % INPUTS: 
        % - typeIn: type of FD
        % - XrefIn: point on which the gradients will be calculated
        % - stepsIn: step(s) used for FD
        % - funIn: handle function (@(x) ...)
        function obj=gradFD(typeIn,XrefIn,stepsIn,funIn)
            %activate or not demo mode
            if nargin > 0
                demo=false;
                if isempty(typeIn)
                    demo=true;
                else
                    obj.type=typeIn;
                end
                if nargin>1
                    obj.Xref=XrefIn;
                end
                if nargin>2
                    obj.stepsDiff=stepsIn;
                end
                %
                if nargin>3
                    obj.fun=funIn;
                    if demo
                        obj.runDemo;
                    end
                end
            else
               obj.displaySchemes; 
            end
            %
            %
        end
        %% setters
        function set.Xref(obj,XX)
            obj.Xref=obj.loadX(XX);
        end
        function set.stepsDiff(obj,steps)
            if all(steps>0)
                obj.stepsDiff=steps;
            else
                fprintf('Bad stepsize - maintain the prévious value(s):\n')
                fprintf('%d',steps);
                fprintf('\n');
            end
        end
        %% getters
        function XX=get.XevalG(obj)
            XX=obj.geneXG();
        end
        function ZZ=get.ZevalG(obj)
            if isempty(obj.fun)&&isempty(obj.ZevalGext)      
                fprintf('Unable to compute the gradients\n Load responses or define a function\n');
            end
            if isempty(obj.fun)&&~isempty(obj.ZevalGext)
                ZZ=obj.ZevalGext;
            end
            if ~isempty(obj.fun)
                ZZ=obj.GcomputeZ;
            end
        end
        function GZ=get.GZeval(obj)
            GZ=obj.computeGZ;
        end
        function stepsOut=get.stepsDiffInternal(obj)
            stepsOut=obj.loadStepsDiff(obj.stepsDiff);
        end
        %% load sample points
        function Xout=loadX(obj,XX)
            Xout=XX;
            obj.dim=size(XX,2);
            obj.nS=size(XX,1);
        end
        %% load the steps for FD
        function stepsOut=loadStepsDiff(obj,steps)
            % deal with specific form of stepsize obj.stepsDiff
            % on row: specific stepsize per dimension
            % on column: specific stepsize per point
            sSteps=size(steps);
            nbR=[1 1];
            if sSteps(1) == 1 
                nbR(1)=obj.nS;
            elseif sSteps(1) ~=1 && sSteps(1) ~= obj.nS
                fprintf(['Wrong size of the steps for computing finite differences (' mfilename ')\n']);
            end
            if sSteps(2) == 1 
                nbR(2)=obj.dim;
            elseif sSteps(2) ~=1 && sSteps(2) ~= obj.dim
                fprintf(['Wrong size of the steps for computing finite differences (' mfilename ')\n']);
            end            
            stepsOut=repmat(steps,nbR);
        end
        %% build coordinates for evaluating the function for gradient
        function XX=geneXG(obj)
            %load configuration
            obj.confFD=loadFD(obj.type);
            %build combination of steps for eveyr dimension
            nbStep=numel(obj.confFD.grad.steps);
            stepsXtmp=[obj.confFD.grad.steps(ones(1,obj.dim),:).';zeros(nbStep*(obj.dim-1),obj.dim)];
            stepsXraw = arrayfun(@(i) circshift(stepsXtmp(:, i), nbStep*(i-1)), 1:obj.dim, 'UniformOutput', false);
            stepsXraw = cell2mat(stepsXraw);
            %build array of points for evaluating the function
            sDiff=obj.stepsDiffInternal;
            XX=zeros(obj.nS*nbStep*obj.dim,obj.dim);
            for itS=1:obj.nS
                nbT=nbStep*obj.dim;
                itX=nbT*(itS-1)+(1:nbT);
                %
                XX(itX,:)=repmat(obj.Xref(itS,:),[nbT 1])+bsxfun(@times,stepsXraw,sDiff(itS,:));
            end
            %remove duplicate and store positions
            [XX,~,obj.dupX]=unique(XX,'rows');
            obj.nX=size(XX,1);
        end
        %% load external Z (external evaluation of the function)
        function loadZextG(obj,ZZ)
            if ~isempty(ZZ)
                if numel(ZZ)==obj.nX
                    obj.ZevalGext=ZZ(:);
                else
                    fprintf('Wrong size of external responses (expected: %i, provided: %i\n',obj.nX,numel(ZZ));
                end
            end
        end
        %% compute responses of the function at the Xeval points
        function ZZ=GcomputeZ(obj,XX)            
            if ~isempty(obj.fun)
                if nargin>1
                    ZZ=feval(obj.fun,XX);
                else
                    ZZ=feval(obj.fun,obj.XevalG);
                end
            else
                fprintf(['Undefined function for evaluation (' mfilename ')\n']);
            end
        end
        %% compute gradients from responses at the Xeval points
        function GZ=computeGZ(obj)
            %build the right Zeval vector
            ZZevalG=obj.ZevalG;
            ddupX=obj.dupX;
            rZeval=repmat(ZZevalG(ddupX),[1,obj.dim]);
            %load coef and divisor
            coefG=obj.confFD.grad.coef;
            nbCoef=numel(coefG);
            divG=obj.confFD.grad.div;
            %build array of coefficients
            coefTmp=[coefG(ones(1,obj.dim),:).';zeros(nbCoef*(obj.dim-1),obj.dim)];
            coefRaw = arrayfun(@(i) circshift(coefTmp(:, i), nbCoef*(i-1)), 1:obj.dim, 'UniformOutput', false);
            coefRaw = cell2mat(coefRaw);
            %product coef*response
            prodZCoef=rZeval.*repmat(coefRaw,[obj.nS,1]);
            %stepsizes
            sDiff=obj.stepsDiffInternal;
            %build the array of gradients
            GZ=zeros(obj.nS,obj.dim);
            for itS=1:obj.nS
                nbT=nbCoef*obj.dim;
                itX=nbT*(itS-1)+(1:nbT);
                GZ(itS,:)=sum(prodZCoef(itX,:),1)./(divG*sDiff(itS,:));
            end
        end
        %% show available schemes
        function displaySchemes(obj)
           %load list of available finite differences
            listT=fieldnames(loadFD); 
            fprintf('List of the available FD technics\n')
            for itT=1:numel(listT)
                fprintf([listT{itT} '\n']);
            end
        end
        %% run demo mode
        function runDemo(obj)
            %load list of available finite differences
            listT=fieldnames(loadFD);
            GZdemo=cell(1,numel(listT));
            for itT=1:numel(listT)
                obj.type=listT{itT};
                GZdemo{itT}=obj.computeGZ;
            end
            fprintf('Display of the results from the available FD technics\n')
            for itT=1:numel(listT)
                fprintf('%s: ',listT{itT});
                if size(GZdemo{1},1)>1;fprintf('\n');end
                for iTG=1:size(GZdemo{itT},1)
                    fprintf('%d ',GZdemo{itT}(iTG,:));
                    fprintf('\n');
                end
            end
        end
    end
    
end
