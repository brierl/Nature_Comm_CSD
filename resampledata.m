function [data2, currHz]=resampledata(data,framerate,tHz,tol)

% resampledata() is an alternative to MATLAB's Signal Processing Toolbox
% resample() command. It has three advantages:
% 
%    1. resample() assumes that the endpoints of your data go to zero.
%    resampledata() assumes that your output should have the same end
%    points as the input. It thus subtracts this linear fit from your data
%    before calling MATLAB's resample().
%    2. resampledata() allows you to change the tolerance with which it
%    resamples.
%    3. resampledata() allows you to input a struct, then resamples each
%    field in the struct. 
% 
% The syntax is:
% 
% [data2]=resampledata(data,framerate,tHz,tol)
% 
% data is your input data, and can be in any format, as long as time is the
% last dimension. framerate is the framrate of the input data. tHz is the
% desired framerate of the output data. tol is the tolerance with which
% resampledata() tries to achieve this exact output framerate. tol is an
% optional input and defaults to 10^{-5}. (If you want to resample a long
% time series to an exact time, use tol around 10^{-5}. If you plan to
% block-average data, then use tol around 10^{-3} so there is a better
% chance that each inter-stimulus interval is resampled to the same
% length.)
% 
% data2 is the same shape as data, with the exception of its length along
% the last dimension. If data is a struct, its fields are individually
% resampled as above and output as fields of data2. 

% Default values
if nargin<4
    tol=1e-5;
end

% This is the way to turn off resampling
if tHz==0
    data2=data;
    currHz=framerate;
    return
else
    currHz=tHz;
end
    
% Check whether you have a struct (then run iteratively) or a matrix (just
% resample directly)
if isstruct(data) % if struct
    % get all fields
    N=fieldnames(data);
    
    % loop through fields
    for f=1:size(N,1)
        fname=N{f};
        % call iteratively
        data2.(fname)=resampledata(data.(fname),framerate,tHz,tol);
    end
else % is a matrix
    % approximate desired resampling ratio as a fraction
    [N,D]=rat(tHz/framerate,tol);
    
    % reshape to everything x time
    [data, Sin, Sout]=datacondition(data,1); 
    Tin=Sout(end); % original length in time

    for n=1:Sout(1) % for every time trace
        signal=squeeze(data(n,:)); % this particular time trace
        
        sigstart=signal(1); % left end point
        sigend=signal(end); % right end point
        
        % linear fit to end points
        alpha1=(sigstart-sigend)/Tin;
        beta=-sigstart;
        
        % remove linear fit
        corrsig=signal+alpha1*(1:Tin)+beta;
        corrsig=double(corrsig);
        % resample with endpoints pinned to zero
        rawresamp=resample(corrsig,N,D);
        
        % re
        Tout=size(rawresamp,2);

        % linear fit to new line
        alpha2=(sigstart-sigend)/Tout;
        
        % add line back in to refit endpoints
        corrresamp=rawresamp-alpha2*(1:Tout)-beta;
        
        data2(n,:)=corrresamp;
    end

    % resampled time
    Tout=size(data2,2);
    
    data2=reshape(data2,[Sin(1:(end-1)) Tout]); % reshape to original shape, except with new time

end

end