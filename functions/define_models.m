function models = define_models(model_name, params_array)
% DEFINE_MODELS Constructs a library of model structs.
%   models = DEFINE_MODELS(model_name, params_array) creates a cell array
%   of model structs for a given model name and an array of parameter structs.
%
%   Inputs:
%   - model_name: A string for the model type (e.g., "Gamma").
%   - params_array: A struct array where each element contains the
%                   parameters for a single model instance.
%                   Example for "Gamma": params(1).alpha = 1; params(1).beta = 10;
%
%   Output:
%   - models: A cell array of final model structs.

    models = {};

    for i = 1:length(params_array)
        params = params_array(i);
        switch model_name
            case "Gamma"
                meantime = params.alpha * params.beta;
                model = struct('name', "Gamma", ...
                    'parameters', params, ...
                    'meantime', meantime, ...
                    'detail_name', sprintf("G(a=%.2g,b=%.2g)", params.alpha, params.beta));
                
            case "Poisson"
                meantime = 1 / params.rate;
                model = struct('name', "Poisson", ...
                    'parameters', params, ...
                    'meantime', meantime, ...
                    'detail_name', sprintf("E(r=%.2g)", params.rate));

            case {"InverseGaussian", "IG", "ig"}
                meantime = params.mu;
                model = struct('name', "InverseGaussian", ...
                    'parameters', params, ...
                    'meantime', meantime, ...
                    'detail_name', sprintf("IG(mu=%.2g,l=%.2g)", params.mu, params.lambda));

            case "deterministic"
                meantime = params.tau;
                model = struct('name', "deterministic", ...
                    'parameters', params, ...
                    'meantime', meantime, ...
                    'detail_name', sprintf("D(t=%.2g)", params.tau));
                
            otherwise
                warning('Model type "%s" not recognized. Skipping.', model_name);
                continue;
        end
        
        models{end+1} = model;
    end
end
