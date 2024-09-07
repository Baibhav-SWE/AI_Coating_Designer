function ref = nload( file, lam )
    %% FUNCTION NLOAD
    % interpolates n&k data from file to a specific wavelength range
    % -------------------------------------------------------------------------
    %% INPUT
    % -------------------------------------------------------------------------
    % file   - 1. name of file with n&k data, three columns: wavelength (nm) | n | k
    %          2. constant real number if starts with "const="
    %          3. n=1 if set to "air"
    % lam    - wavelength (nm) for interpolation
    % -------------------------------------------------------------------------
    %% OUTPUT
    % -------------------------------------------------------------------------
    % ref - interpolated complex refractive index: ref = n + ik
    % -------------------------------------------------------------------------
    switch file
        case "air"
            ref = ones(1,numel(lam));
        otherwise 
            if startsWith(file,"const=") && ~contains(file,"j")
                ref = str2double(extractAfter(file,"const="))*ones(1,numel(lam));
            elseif startsWith(file,"const=") && contains(file,"j")
                realn = str2double(extractBetween(file,"const=","+"));
                imagn = str2double(extractBetween(file,"+","j"));
                ref = complex(realn,imagn)*ones(1,numel(lam));
            else
                data = load(file);
                v = fieldnames(data);
                ntable = data.(v{1});
                ref = complex( interp1( ntable(:,1), ntable(:,2), lam, 'spline' ),...
                             interp1( ntable(:,1), ntable(:,3), lam, 'spline') );
            end
    end
    % -------------------------------------------------------------------------
    end
    % -------------------------------------------------------------------------