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
        XevalG;             % points used to compute FD
        ZevalG;             % responses used to compute FD
        GZeval;             % values of gradients
        XevalH;             % points used to compute hessians with FD
        ZevalH;             % responses used to compute hessians with FD
        HZeval;             % values of hessians
        fun;                % considered function
        stepSizes=1e-4;     % steps of FD
        nbV;                % dimension (number of design variables)
        nS;                 % number of sample points
        nX;                 % number of points used to compute FD (duplicate points are removed)
        lb;                 % lower bound
        ub;                 % upper bound
        partialData;        % data if partial gradients calculation
    end
    properties (Access=private)
        confFD;             % configuration of the chosen scheme used for FD
        dupX;               % duplicate coordinates of Xeval (for reducing the number of evaluations of the function
        stepsDiffInternal;  % steps of FD
        ZevalGext;          % array of responses used in the case of external function
        demoMode=false;     % flag for demo mode
        buildOk=false;      % flag to check if data has been built
        maskParaVal={};     % mask to existing current sample point
        valParaGrad={};     % values of sample points
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% constructor of finite differences class
        function obj=gradFD(varargin)
            %activate or not demo mode
            if nargin > 0
                obj.loadClass(varargin{:})
            else
                obj.displaySchemes;
            end
            %
            %
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% setters
        function set.Xref(obj,XX)
            obj.Xref=obj.loadX(XX);
        end
        function set.type(obj,valI)
            if size(valI,2)<3
                valI(:,end+1:3)=' ';
            elseif size(valI,2)>3
                valI(:,3:end)=[];
            end
            %
            obj.type=valI;
            obj.resetFlag;
        end
        function set.stepSizes(obj,valI)
            if all(valI>0)
                if numel(valI)==1
                    valI=repmat(valI,1,obj.nbV);
                end
                obj.stepSizes=valI;
                obj.resetFlag;
            else
                fprintf('Bad stepsize - maintain the previous value(s):\n')
                fprintf('%d',obj.stepSizes);
                fprintf('\n');
            end
        end
        function set.lb(obj,valI)
            obj.lb=valI;
            obj.resetFlag;
        end
        function set.ub(obj,valI)
            obj.ub=valI;
            obj.resetFlag;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % reset build flag
        function resetFlag(obj)
            obj.buildOk=false;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % build if necessary
        function buildCond(obj)
            if ~obj.buildOk
                obj.buildData;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% getters
        function XX=get.XevalG(obj)
            XX=obj.points();
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
            stepsOut=obj.loadStepsDiff(obj.stepSizes);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % method to load the class
        % INPUTS:
        % - typeIn: type of FD
        % - XrefIn: point on which the gradients will be calculated
        % - stepsIn: step(s) used for FD
        % - funIn: handle function (@(x) ...)
        % - lbIn, ubIn: lower/upper bounds
        % - partialIn: specify which gradients have to be calculated 
        function loadClass(obj,typeGradIn,paraValIn,stepSizesIn,funIn,lbIn,ubIn,partialIn)
            
            %number of variables
            obj.nbV=numel(paraValIn);
            %deal with input arguments
            if nargin>1
                obj.type=typeGradIn;
            end
            if nargin>2
                obj.Xref=obj.loadX(paraValIn);
            end
            if nargin>3
                obj.stepSizes=stepSizesIn;
            end            
            %default values
            lbDef=-Inf*ones(obj.nbV,1);
            ubDef=+Inf*ones(obj.nbV,1);            
            %
            if nargin>5
                obj.lb=lbDef;
                obj.ub=ubDef;
                if ~isempty(lbIn)
                    obj.lb=lbIn;
                end
            end
            if nargin>6            
                if ~isempty(ubIn)
                    obj.ub=ubIn;
                end
            end
            %
            if nargin>7
                obj.partialData=partialIn;
            end
            %
            if nargin>4
                obj.fun=funIn;
                if obj.demoMode
                    obj.runDemo;
                end
            end
            %
            obj.buildOk=false;            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% load sample points
        function Xout=loadX(obj,XX)
            Xout=XX;
            obj.nbV=size(XX,2);
            obj.nS=size(XX,1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% load the steps for FD
        function stepsOut=loadStepsDiff(obj,steps)
            % deal with specific form of stepsize obj.stepSizes
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
                nbR(2)=obj.nbV;
            elseif sSteps(2) ~=1 && sSteps(2) ~= obj.nbV
                fprintf(['Wrong size of the steps for computing finite differences (' mfilename ')\n']);
            end
            stepsOut=repmat(steps,nbR);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % adapt type per points and dimension
        function adaptType(obj)
            %scheme per sample point and dimension
            if size(obj.type,1)==1           
                typeRaw=obj.type;
                obj.type=cell(obj.nS,obj.nbV);
                for itS=1:obj.nS
                    for itV=1:obj.nbV
                        obj.type{itS,itV}=typeRaw;
                    end
                end
            end    
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % build data for computing gradients
        function buildData(obj)            
            %depending on the situation (bounded or partial gradients
            %computation)
            if ~isempty(obj.partialData)||any(~isinf(obj.ub))||any(~isinf(obj.lb))
                %adapt  type of FD scheme
                obj.adaptType;
                obj.buildDataConstrained;
            else
                obj.buildDataUnconstrained;
            end               
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % build data for computing gradients with constraint(s)
        function buildDataPerPoint(obj)
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % build data for computing gradients with constraint(s)
        function buildDataConstrained(obj)
            %along the sample points
            for itS=1:obj.nS
               %obj.buildDataPerPoint(obj.Xref(itS,:),) 
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % build data for computing gradients with no constraint
        function buildDataUnconstrained(obj)
        %load configuration
            obj.confFD=loadFD(obj.type);
            %build combination of steps for every dimension
            nbStep=numel(obj.confFD.grad.steps);
            stepsXtmp=[obj.confFD.grad.steps(ones(1,obj.nbV),:).';zeros(nbStep*(obj.nbV-1),obj.nbV)];
            stepsXraw = arrayfun(@(i) circshift(stepsXtmp(:, i), nbStep*(i-1)), 1:obj.nbV, 'UniformOutput', false);
            stepsXraw = cell2mat(stepsXraw);
            %build array of points for evaluating the function
            sDiff=obj.stepsDiffInternal;
            XX=zeros(obj.nS*nbStep*obj.nbV,obj.nbV);
            for itS=1:obj.nS
                nbT=nbStep*obj.nbV;
                itX=nbT*(itS-1)+(1:nbT);
                %
                XX(itX,:)=repmat(obj.Xref(itS,:),[nbT 1])+...
                    bsxfun(@times,stepsXraw,sDiff(itS,:));
            end
            %remove duplicate and store positions
            [obj.valParaGrad,~,obj.dupX]=unique(XX,'rows');
            obj.nX=size(obj.valParaGrad,1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % export points to calculate responses
        function ptsVal=points(obj)
            %
            obj.buildCond;
            %
            if iscell(obj.valParaGrad)
                paraValV=vertcat(obj.valParaGrad{:});
                ptsVal=paraValV(~vertcat(obj.maskParaVal{:}),:);
            else
                ptsVal=obj.valParaGrad;
            end
        end
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% compute gradients from responses at the Xeval points
        function GZ=computeGZ(obj)
            %build the right Zeval vector
            ZZevalG=obj.ZevalG;
            ddupX=obj.dupX;
            rZeval=repmat(ZZevalG(ddupX),[1,obj.nbV]);
            %load coef and divisor
            coefG=obj.confFD.grad.coef;
            nbCoef=numel(coefG);
            divG=obj.confFD.grad.div;
            %build array of coefficients
            coefTmp=[coefG(ones(1,obj.nbV),:).';zeros(nbCoef*(obj.nbV-1),obj.nbV)];
            coefRaw = arrayfun(@(i) circshift(coefTmp(:, i), nbCoef*(i-1)), 1:obj.nbV, 'UniformOutput', false);
            coefRaw = cell2mat(coefRaw);
            %product coef*response
            prodZCoef=rZeval.*repmat(coefRaw,[obj.nS,1]);
            %stepsizes
            sDiff=obj.stepsDiffInternal;
            %build the array of gradients
            GZ=zeros(obj.nS,obj.nbV);
            for itS=1:obj.nS
                nbT=nbCoef*obj.nbV;
                itX=nbT*(itS-1)+(1:nbT);
                GZ(itS,:)=sum(prodZCoef(itX,:),1)./(divG*sDiff(itS,:));
            end
        end
        function GZ=computeGrad(obj,respIn,currentResp)
            if nargin>1
                obj.loadZextG(respVal);
                %
                if nargin>2
                    
                end
            end
            obj.loadZextG(respVal);
            %
            GZ=obj.computeGZ;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% show available schemes
        function displaySchemes(obj)
            %load list of available finite differences
            listT=fieldnames(loadFD);
            fprintf('List of the available FD schemes\n')
            for itT=1:numel(listT)
                fprintf([listT{itT} '\n']);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
