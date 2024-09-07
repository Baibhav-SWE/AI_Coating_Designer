function stack = set_stack( mat, thick, lam, angle, varargin )
    %% FUNCTION set_stack
    % creates dictionary for the multilayered stack
    % -------------------------------------------------------------------------
    %% INPUT
    % -------------------------------------------------------------------------
    % mat   - string array with names for files containing n&k for materials
    % thick - thicknesses of layers, notice first and last layer are infinite
    % lam   - wavelengths
    % angle - angles of inicdence
    % varargin - optional arguments:
    %            1. incoh - threshold thickness for incoherence in layers
    % -------------------------------------------------------------------------
    %% OUTPUT
    % -------------------------------------------------------------------------
    % stack - dictionary of stack ready for further analysis
    % -------------------------------------------------------------------------
    %% MODIFYING INPUT IF NEEDED
    % -------------------------------------------------------------------------
    if size(lam,1) == 1
        lam = lam';
    end
    % -------------------------------------------------------------------------
    %% INTERPOLATING N&K VALUES
    % -------------------------------------------------------------------------
    [mat_unique,~,ic] = unique(mat);
    n_unique = zeros(numel(lam),numel(mat_unique));
    %
    for imat = 1:numel(mat_unique)
        n_unique(:,imat) = nload(mat_unique(imat),lam);
    end
    n = n_unique(:,ic);
    % -------------------------------------------------------------------------
    %% CREATING STACK DICTIONARY
    % -------------------------------------------------------------------------
    if ( length(varargin) == 1 )
        nincoh = thick<varargin{:};
    else 
        nincoh = true(size(thick));
    end
    keys = ["nk","thick","nincoh","wavelength","angle"];
    values = {n, thick, nincoh, lam, angle};
    % -------------------------------------------------------------------------
    stack = dictionary(keys,values);
    % -------------------------------------------------------------------------
    end
    % -------------------------------------------------------------------------