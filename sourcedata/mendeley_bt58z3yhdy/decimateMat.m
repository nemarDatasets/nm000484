function odata = decimateMat( idata, r )
%%
%  This function follows the method of matlab's decimate function, but
%  here it can also handle non-integer 'r' values.  
% 

  assert( size(idata,1) < size(idata,2) );
    nd = size(idata,2);
    m = size(idata,1);
    nout = ceil(nd/r);
    nfilt = 8;
    rip = .05;	% passband ripple in dB
    [b,a] = cheby1(nfilt, rip, .8/r);
    while all(b==0) || (abs(filtmag_db(b,a,.8/r)+rip)>1e-6)
        nfilt = nfilt - 1;
        if nfilt == 0
            break
        end
        [b,a] = cheby1(nfilt, rip, .8/r);
    end
    if nfilt == 0
        error(message('signal:decimate:InvalidRange'))
    end
    % be sure to filter in both directions to make sure the filtered data has zero phase
    % make a data vector properly pre- and ap- pended to filter forwards and back
    % so end effects can be obliterated.
    odata = filtfilt(b,a,idata')';
    nbeg = floor(r - (r*nout - nd));

    if( r == floor(r) )
        odata = odata(:,nbeg:r:nd);  % integer value after all...
    else
        idx = [nbeg:r:nd];
        lower = floor(idx);
        upper = ceil(idx);
        w     = idx - lower;

        if( floor(nbeg) == 0 )
          lower(lower<1) = 1;
          upper(upper<1) = 1;
        end
        n_odata = size(odata,2);
        if( ceil(nd) > n_odata )
          lower(lower>n_odata) = n_odata;
          upper(upper>n_odata) = n_odata;
        end
        odata = odata(:,lower).*repmat((1-w),m,1) + odata(:,upper).*repmat(w,m,1);
    end
    if( diff(size(idata))*diff(size(odata)) < 0 )
        odata = odata';
    end
end

%--------------------------------------------------------------------------
function H = filtmag_db(b,a,f)
%FILTMAG_DB Find filter's magnitude response in decibels at given frequency.

    nb = length(b);
    na = length(a);
    top = exp(-1i*(0:nb-1)*pi*f)*b(:);
    bot = exp(-1i*(0:na-1)*pi*f)*a(:);
    
    H = 20*log10(abs(top/bot));
end
