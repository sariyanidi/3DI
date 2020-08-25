% The images will by resized by resize_coef prior to processing. To speed
% up the code, this can be reduced up to ~0.25 but results will be slightly
% worse
resize_coef = 0.65;
files = dir('./input/AFLW3D/*png');


for fi=1:length(files)
    bcoef = 1;
    
    fpath = sprintf('./input/AFLW3D/%s', files(fi).name);
    [~, bn] = fileparts(fpath);
    fprintf('%03d) Processing %s... ', fi, bn);
    
    % We may need to attempt fitting multiple times:
    % if the inequality constraints for landmarks are not satisfied, we
    % will slightly enlarge the interval of the landmark locations (with bcoef) 
    % and will try again.
    %
    % There are multiple reasons to why constraints may not be satisfied at
    % first try. One reason is that the we assume lens distortion to be 
    % negligible, and if it's not, the warping that it causes may not be 
    % accurately reconstructed by the 3D morphable model. 
    for attempt=1:5
        clear result
        
        try
            tic
            % The code that actually does the fitting is fit3dmm_to_image()
            % 
            % The fitting parameters are stored in fit_params
            % Check the comments in the render_3dhead file to see how 
            % the fitting parameters can be used to reconstruct 
            % the dense 2D or 3D face shape as well as the facial landmarks
            [data, opts, fit_params, result, existed] = fit3dmm_to_image(fpath, bcoef, resize_coef);
            
            if ~existed
                fprintf('Took %.2f secs\n', toc);
                render_3dhead(data, opts, fit_params);
            else
                dev_null = toc;
            end
            
            if result
                break
            end
        catch e 
            fprintf('%s\n', e.message);
            bcoef = bcoef*1.2;
        end
    end
end
