function features = compute_feature_demo( iEEGfileName, hfoFileName, fs, deciFs,filtRange,secWin,secBak,secBUFFER)
  %
  %  Compute wavelet entropy for qHFO detections in a given file
  %
  %

  %% compute number of samples 
  nWin = floor(secWin*deciFs); % length of the waveform that is used to construct time-frequency map via wavelet transformation
  nBak = floor(secBak.*deciFs); % length of the waveform that is used estimate wavelet entropy
  BUFFER = round(secBUFFER*fs); % length of the buffer used to mitigate ringing artifacts near the edges of the windowed signals

  %% construct filter
  if( ~isnumeric(filtRange))
      assert(strcmp(filtRange,'raw'),'only expect ''raw''as str input');
      doFilt = false;
  else
      assert(numel(filtRange) == 2);
      doFilt = true;
  end
  if doFilt
      Hd = [];
      if filtRange(1)>0
          Hd = cauerfilt( 80, min( 500, 0.4*deciFs ), deciFs );
      else
          assert(filtRange(1)==0,'unexpected filter setting');
          [z, p, k ] = butter( 10, filtRange(2) / (deciFs/2) );
          [Hd.sosMatrix, Hd.ScaleValues] = zp2sos( z, p, k );
      end
  end

  %% first timer
  ticA = tic;
  times.load = 0;
  times.process = 0;
  
  %% load data
  ieeg = matfile(iEEGfileName); 
  nSamples = size(ieeg.data,2);
  hfo = load(hfoFileName);

  %% initialize the vector to store features
  features = nan(length(hfo.startIdx),1);

  %% compute locations(k1,k2) of HFO detections with buffer
  buf_1 = ones(size(hfo.startIdx)) * BUFFER;
  I = hfo.startIdx <= BUFFER;
  buf_1(I) = hfo.startIdx(I) - 1;

  buf_2 = ones(size(hfo.stopIdx)) * BUFFER;
  I = hfo.stopIdx >= nSamples - BUFFER;
  buf_2(I) = nSamples - hfo.stopIdx(I);
  
  k1 = double(hfo.startIdx) - buf_1;
  k2 = double(hfo.stopIdx) + buf_2;
  
  I = k2 - k1 + 1 > 2*BUFFER + 5 & k2 - k1 > 60; % skip HFOs that  are too short for filtering

  %% loop
  for i=find(I(:)')
      % load
      ticB = tic;
      x = ieeg.data(:, k1(i): k2(i) );
      times.load = times.load + toc(ticB);
      ticB = tic;
      x = double(x(hfo.chanIdx(i),:)) - mean(x); % common average refererence
      filtOK = true;
      %% filter if needed
      if doFilt
          filtOK = false;
          try
              x = decimateMat( x, fs/deciFs );
              x = filtfilt( Hd.sosMatrix, Hd.ScaleValues, x );
              filtOK = true;
          catch err
              display(err)
              display(err.message)
              display(err.stack(1))
          end
      end
      if( filtOK )
        % trim
        buf_1(i) = round(buf_1(i)/fs*deciFs); 
        buf_2(i) = round(buf_2(i)/fs*deciFs); 
        hfoL =  length(x)-buf_1(i)-buf_2(i); %% length of hfo detection
   
	if hfoL < nWin      
            midI = floor(hfoL/2)+buf_1(i); %% center of HFO detection
            I = [ (midI-ceil(nWin/2)+1):( midI+floor(nWin/2) ) ]; %% index of the signal that is used to construct time-frequency mapvia wavelet transformation. This signal is centered around the detected HFO.
            bakI = [ (midI-ceil(nBak/2)+1):(midI+floor(nBak/2)) ]-I(1)+1; %% index of the signal of interest that is centered around the detected HFO
	else
            bakI = 1:nBak;
	end
        try
          features(i,:) = single_waveform_feature_demo(x(I),bakI, deciFs); % compute wavelet entropy for current HFO detection
        catch err
            display(err)
            display(err.message)
            display(err.stack(1))
        end

      end
  
      times.process = times.process + toc(ticB);
      
  end
  
  ticB = tic;
  times.total = toc(ticA);
  
  fprintf('Processing times:\n');
  fprintf('\tload:    %7.3f\n', times.load );
  fprintf('\tprocess: %7.3f\n', times.process );
  fprintf('\ttotal:   %7.3f\n', times.total );
  

