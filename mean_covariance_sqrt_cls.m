classdef mean_covariance_sqrt_cls
    %A container for the mean and a square root of covariance matrix
    properties
        dim %dimension
        data%of n*(n+1) dimension, stores mean and a squareroot
    end
    methods
        function obj = mean_covariance_sqrt_cls(x_or_xS,Sigma_sqrt)
            %overloaded constructor: either input xS in a single matrix, or
            %separately by its center and a sqrt of covariance matrix.
            if nargin == 1
                obj.data = x_or_xS;
                obj.dim = size(x_or_xS,1);
            elseif nargin == 2
                obj.dim = numel(x_or_xS);
                obj.data(:,1) = x_or_xS(:);
                obj.data(:,2:obj.dim+1) = Sigma_sqrt;
            else
                error('wrong number of input arguments')
            end
        end
        function [mu] = mean(obj)
            mu = obj.data(:,1);
        end
        function [Sigma] = covariance(obj)
            M = obj.data(:,2:obj.dim+1);
            Sigma = M*M';
        end
        function [M] = c_sqrt(obj)
            M = obj.data(:,2:obj.dim+1);
        end
        function plot_ellipse(obj,proj_dimension,varargin)
            %Plots covariance ellipse projected into specified dimensions. 
            % Passes varargin to plot function
            mu = obj.mean();
            Sigma = obj.covariance();
            mu_p = mu(proj_dimension);
            Sigma_p = Sigma(proj_dimension);
            N = 100;
            t = linspace(0,2*pi,N);
            [V,D] = eig(Sigma_p);
            U0Orth = sqrt(diag(D))'.*V;
            xy = mu_p + U0Orth(:,1)*cos(t) + U0Orth(:,2)*sin(t);
            if nargin==1
                plot(xy(1,:),xy(2,:),'k')
            else
                plot(xy(1,:),xy(2,:),varargin{:});
            end
        end
    end
end