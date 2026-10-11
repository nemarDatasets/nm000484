classdef EMGdetClass < handle

    properties( SetAccess = protected )

      
      %%input created in the load_par function 
      par
      %% parameters that will be from par.XXX
      fs % sampling rate of the data feed into EMG detector
      thresAcp  % Acp refers to the LineLength threshold that will be "accepted" as within the physiological range. Below = physiological
      thresDet  % Det is the threshold over which a signal is detected as EMG.  Above = suspected EMG
      thresLR   % ratio between line length and range of amplitude.   Above = EMG candidate
      agreeNum  % checks whether at least Num channels have simultaneous detections
        
        
      %% parameters that depend on input data
      FS %% sampling rate in raw data
      buffer %% length of buffer used to eliminate ringing effect of filtering around the start/stop of signals
      chanNum % number of channels
      epSam  %number of samples within an epoch of 'epoch' (in function DoAll) minutes
      wdSam  %number of samples within a window of length par.wdSec seconds
      shSam  %number of samples of the sliding window step par.shSec
      blSam  %number of samples within the baseline of par.bsSec seconds
      blWdNum  %number of windows within the baseline
      shiftPerWd % number of steps within a window
      epWdNum %number of windows within an epoch
      epEffSam %number of effective samples within an epoch after truncating off extra samples shorter than a window at the end of the epoch
   
      wdPerShift  % minimum number of windows that can cover the entire effective epoch if there's no overlapping
      tailWd % number of extra steps exceed the epoch according to wdPerShift
      shEpSam % number of samples in a temporary vector used to compute the average line-length values from overlaping windows 
      resetBaseline % whether or not to reset the baseline when start processing a new epoch

      %%function handles
      compValue % compute relative line-length (for detemining EMG artifect) in comparison with baseline 
      compValueR % compute range
      compValueL % compute line-length
      compAcp % compute relative line-length (for detemining physiological variation) in comparison with baseline

      %%comb filter
      haFilt
      
      %%data 
      X %% EEG signals (channel by sample) 
      baseline %% baseline line length values
      lenSam %% line length values indexed by sample
      ranSam %% range values indexed by sample
      markerSam %% binary EMG markings indexed by sample
      startIdx %% start index of EMG detection
      stopIdx %% stop index of EMG detection
      chanIdx %% channel index of EMG detection
      startIdxDe  %% start index of EMG detection with respect to decimated data
      stopIdxDe %% stop index of EMG detection with respect to decimated data

      %% flags indicating status of prepossing steps
      isDecimate %% decimated ?
      isFilt %% filtered ?
      comAve %%  common average is applied ?
    end
    methods
      function obj = EMGdetClass( par )
	
	      %parameters from par  
              obj.fs = par.fs;
              obj.resetBaseline = par.resetBaseline;
	      obj.thresDet = par.thresDet;
              obj.thresAcp = par.thresAcp;
              obj.thresLR = par.thresLR;
              obj.agreeNum = par.agreeNum;

	      %calculate sample numbers
              obj.wdSam = par.wdSec * par.fs ;
              obj.shSam = par.shSec * par.fs ;
              obj.blSam = par.blSec * par.fs ;
              obj.shiftPerWd = floor( obj.wdSam/obj.shSam );
              obj.blWdNum = floor((obj.blSam-obj.wdSam)/obj.shSam)+1;

            % comb filter
              if obj.fs/60 == floor(obj.fs/60)
                d  = fdesign.comb('notch','N,BW',obj.fs/60,0.02); %use matlab's comb filter at 60*N Hz
                Hd = design(d,'SystemObject',true);
                obj.haFilt.b = Hd.Numerator; obj.haFilt.a = Hd.Denominator;
              else
                r=0.95;f0=60;
                [~,obj.haFilt.b,obj.haFilt.a]=harmFilt( r, f0, obj.fs, 0 );
              end

            % functions
              obj.compValueR=@(x)(range(x));%range
              obj.compValueL=@(x)(sum(abs(diff(x))));%lineLength

              switch par.featureName  % decide which method to compute the relative line-length value(for determining EMG artifect) in comparison with baseline.  We used 'lenMu'.
                case 'lenStd'
                  obj.compValue=@(Len,baseline)((Len-mean(baseline))/std(baseline));
                case 'lenMu'
                  obj.compValue=@(Len,baseline)(Len/mean(baseline));
                case 'lenValue'
                  obj.compValue=@(Len,baseline)(Len);
              end
              switch par.acpMethod % decide which method to compute the relative line-length value(for determining physiological variation) in comparison with baseline.  We used 'mu'
                case 'std'
                  obj.compAcp=@(Len,baseline)((Len-mean(baseline))/std(baseline));
                case 'mu'
                  obj.compAcp=@(Len,baseline)(Len/mean(baseline));
              end

        end

        function doAll(obj, data,FS, baseline,epoch,buffer)
	  %% main function to perform EMG detection
          obj.initialize(data,FS,baseline,epoch,buffer); % initialize variables based on input information 
          obj.prepData(); % preprocessing
          obj.getValue(); % compute relative line length and range values
          obj.getDet(); % get EMG detections
        end

        function initialize(obj,data,FS,baseline,epochSec,buffer)

	  %% initalize variables from input data
          obj.buffer = buffer;
          obj.FS = FS;
          obj.X = data;

	  %% compute sample numbers based on length of input data
          obj.epSam = floor( epochSec * obj.fs );
          obj.epWdNum = floor((obj.epSam-obj.wdSam)/obj.shSam)+1;
          obj.epEffSam = floor( (obj.epSam-obj.wdSam)/obj.shSam )*obj.shSam+obj.wdSam;
          obj.wdPerShift = ceil(obj.epWdNum/obj.shiftPerWd);
          obj.tailWd = obj.wdPerShift * obj.shiftPerWd-obj.epWdNum;
          obj.shEpSam = (obj.wdPerShift+1)*obj.wdSam-obj.shSam;

	  %% initialize variables indexed by sample
	  obj.lenSam = zeros(obj.chanNum,obj.epEffSam);
          obj.ranSam = zeros(obj.chanNum,obj.epEffSam);
          obj.markerSam = false(obj.chanNum,obj.epEffSam);

	  %% making sure EEG signals on individual channels correspond to rows 
          if diff ( size( obj.X ) ) < 0
            obj.X = obj.X';
          end
          obj.chanNum = size(obj.X,1);
          if obj.resetBaseline
              obj.baseline = [];
          else
              obj.baseline = baseline;
          end

	  %% flags indicating prepossing status
          obj.isFilt = false;
          obj.isDecimate = false;
          obj.comAve = false;
        end

        function getDet(obj)
	    %% identify samples whose line-length and range values meet the criteria for EMG artifact

            obj.markerSam = (obj.lenSam>obj.thresDet) & (obj.ranSam>obj.thresLR); %% EMG candidates
            markerAgr = sum( obj.markerSam ) > obj.agreeNum ; %% indexes of samples with at least N(=obj.agreeNum) simultaneous EMG candidates 
            fprintf('total EMG duration: %0.1f sec, ', sum(markerAgr)/obj.fs);
            obj.markerSam = obj.markerSam &  markerAgr; %% final EMG markings on different channels
            idxTemp = cell(1,obj.chanNum); %% temporary container for start/stop indexes of continuous EMG artifects 
            for c = 1: obj.chanNum
              regionSam = obj.marker2region(obj.markerSam(c,:)); %% extract start/stop index from binary EMG markings  
              idxTemp{c} = [regionSam ; repmat(c,1,size(regionSam,2))]; %% adding channel indexes(3rd row) to start/stop indexes(1st/2nd rows) 
            end
            idxTemp = cell2mat(idxTemp); %% 1st row: start index, 2nd row: stop index, 3rd row: stop index

            fprintf('accumulated EMG duration %0.1f == %0.1f sec\n',sum(sum(obj.markerSam))/obj.fs,...
                    sum( diff(idxTemp(1:2,:)) )/obj.fs );

            obj.startIdxDe = idxTemp(1,:); 
            obj.stopIdxDe = idxTemp(2,:);
            obj.chanIdx = idxTemp(3,:);
            obj.startIdx = floor( idxTemp(1,:)*(obj.FS/obj.fs) );
            obj.stopIdx = ceil( idxTemp(2,:)*(obj.FS/obj.fs) );

        end

        function regionSam = marker2region(obj,markerSam)
	  %% extract start/stop indexes of continuous EMG artifacts from binary EMG markings  
          markerDiff=[0 diff(markerSam)]; %% EMG starts/stops the points where the value of binary markings changes  
          bgIdx=find(markerDiff==1);
          edIdx=find(markerDiff==-1);
          if markerSam(1)==1
            bgIdx=[1 bgIdx];
          end
          if markerSam(end)==1
            edIdx=[edIdx length(markerSam)];
          end
          assert(length(bgIdx)==length(edIdx));
          regionSam = [bgIdx;edIdx];
        end

        function prepData(obj)
	  %% preprocessing steps
	  %% apply common average reference
          if ~obj.comAve
            obj.X = double(obj.X); 
            obj.X = obj.X - mean(obj.X);
            obj.comAve = true;
          end
	  %% decimate data
          if ~obj.isDecimate 
            obj.X = decimateMat( obj.X, obj.FS/obj.fs );
            obj.isDecimate = true;
          end
	  %% apply comb filter to remove 60Hz harmonics
          if ~obj.isFilt
            obj.X = filtfilt( obj.haFilt.b, obj.haFilt.a, obj.X' )';
            obj.isFilt = true;
          end
          buffer = floor ( ( size(obj.X,2) - obj.epSam )/2 );
          obj.X = obj.X(:, buffer + 1 : buffer + obj.epSam );
          %%fprintf('actual buffer length %d, effective buffer length %d\n',obj.buffer,buffer/obj.fs*obj.FS);
          assert( size(obj.X,2)==obj.epSam );  
        end

        function getBaseline(obj)
	  %% initialize baseline values by randomly select line length values from individual channels  
          disp('recomputing baseline...');
          assert ( obj.isDecimate & obj.isFilt );
          obj.baseline = zeros(obj.chanNum,obj.blWdNum);
          N = floor(obj.epSam/obj.wdSam);
          for c = 1:obj.chanNum % iterate over channels
            X=obj.X(c,1:N*obj.wdSam);
            wdX=reshape(X,obj.wdSam,N);
            lenValue=obj.compValueL(wdX);
            obj.baseline(c,:)=randsample( lenValue, obj.blWdNum );
          end
        end

        function getValue(obj)
	  %% compute line length and range values using a sliding window, and then assign those values to data samples
          assert ( obj.isDecimate & obj.isFilt );
          if isempty(obj.baseline)
            obj.getBaseline;
          end
          for c = 1: obj.chanNum %% iterate over all channels
            w = [ 1, obj.wdSam ]; %% location of sliding window
            lenWd = zeros(1,obj.epWdNum); %% line length values of windows 
            ranWd = zeros(1,obj.epWdNum); %% range values of windows
            for i = 1 : obj.epWdNum  %% iterate over sliding windows
              newL = obj.compValueL(obj.X(c,w(1):w(2))); %% line length value of the current window
              newR = obj.compValueR(obj.X(c,w(1):w(2))); %% range value of the current window
              lenWd(i) = obj.compValue( newL, obj.baseline(c,:)); %% relative line length of the current window in comparsion with baseline used for identifying EMG contamination
              ranWd(i) = ( newL - mean(obj.baseline(c,:)) )/newR; %% ratio between line length and range
              acpValue=obj.compAcp(newL,obj.baseline(c,:));  %% relative line length of the current window in comparsion with baseline used for identifying physiological variation
              if acpValue <= obj.thresAcp %% if the relative line length is small enough then the current value is used to update baseline 
		obj.baseline(c,:) = [obj.baseline(c,2:end) newL];
              end
              w=w+obj.shSam; %% move the window one step forward
            end  
            obj.lenSam(c,:) = obj.getValueSam(lenWd); %% assigning relative line length value to corresponding samples
            obj.ranSam(c,:) = obj.getValueSam(ranWd); %% assigning range value to corresponding samples
          end
        end
        function valueSam = getValueSam(obj,valueWd)
	  %% transforming window-indexed values to sample-indexed values, where the final value for each sample is the average of the values from overlaping windows
          valueShift = zeros(obj.shiftPerWd,obj.shEpSam);
          valueWd = [valueWd zeros(1,obj.tailWd)]; %% window-indexed values with zeros tappings at the end
          valueWd = reshape(valueWd,obj.shiftPerWd,[]); %% transform window-indexed values so that overlapping windows tend to be in the adjecent columns 

          for i=1:obj.shiftPerWd
            valueShift(i,:)=[zeros(1,obj.shSam*(i-1)) obj.wd2sample(valueWd(i,:)) zeros(1,obj.shSam*(obj.shiftPerWd-i))]; %% shifitting window-indexed values so that values correpsongding the same samples are in the same column  
          end

          valueSam = mean(valueShift(:,1:obj.epEffSam),1); %% final value for each sample is the average of the values from overlaping windows

        end

        function sVal = wd2sample(obj,wVal)
	  %% transforming window-indexed values to sample-indexed values via replication 
          if size(wVal,1)>1
            wVal=wVal';
          end
            sVal=repmat(wVal,obj.wdSam,1);
            sVal=reshape(sVal,[],1);
            sVal=sVal';
        end

        function Delete(obj)
            delete(obj);
        end

    end%end of methods Def

end%end of class def
