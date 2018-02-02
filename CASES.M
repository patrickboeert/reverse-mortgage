classdef cases
% CASES  dataset of HECM cases

   properties
       
       data = [];   % dataset object
           
   end % properties

   methods
       
      %% CONSTRUCTOR
      function obj = cases()
        % IMPORT  creates dataset of HECM cases from HECM data
        %   
        % creates MATLAB dataset object from HECM case data exported
        % from STATA; it also computes and adds variables relevant for
        % simulation of the portfolio
        %
        % N = 390626           
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Dataset Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %           access: {1,1} for individual element
        %                   (1,1) for observation
        %                   dataset.(varname) for dynamic indexing
        %
        %           date conversion functions:    
        %                    datenum(), datevec(), datestr()
        %           
        %           dataset functions: 
        %                    nominal(), ordinal(), dataset()
        %
        %           nominal class:
        %                    getlabels(), setlabels(), islevel(), levelcounts(),
        %                    ismember(), intersect(), union(), mergelevels(),
        %                    setdiff(), setxor(), unique(), summary(),
        %                    isundefined(),
        % 
        %           also see: help categorical/methodname
        %                     help nominal/methodname
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Unique Codes: Codierungswert (# of Obs)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % For all date variables, NaNs are coded as 00-Jan-0000 in
        % string representation (resp. 0 in elapsed MATLAB time 
        % representation)
        %
        %  1 actn_cd:               1, 2, 15, 21, 23, 30, 60, 70, 71, 74, 99, NaN (339699)
        %  2 actn_dt:               date (elapsed), NaN (339692)
        %  3 adjstmnt_dt:           date (elapsed), NaN (381)
        %  4 borr_dt_of_brth:       date (elapsed), NaN (129435)
        %  5 borr_dt_of_dth:        date (elapsed), NaN (388061)
        %  6 borr_gender:           1, 3, NaN (1892)
        %  7 borrower_age:          double, NaN (261770)
        %  8 clsng_dt:              date (elapsed), NaN (310761)
        %  9 coborr_dt_of_brth:     date (elapsed), NaN (296082)
        %  10 coborr_dt_of_dth:     date (elapsed), NaN (389247)
        %  11 coborr_gender:        1, 3, NaN (250579)
        %  12 cs_estblshd_dt:       date (elapsed), NaN (0)
        %  13 cs_sts:               1, 2, 4, 5, 6, 8, 9, NaN (0)
        %  14 db_dt:                date (elapsed), NaN (1)
        %  15 endrsmnt_dt:          date (elapsed), NaN (0)
        %  16 flr_rt:               double, NaN (0)
        %  17 fndng_dt:             date (elapsed), NaN (1)
        %  18 init_fee:             double, NaN (0)
        %  19 init_prncpl_lmt:      double, NaN (0)
        %  20 int_rt:               double, NaN (0)
        %  21 int_rt_10yr:          double, NaN (0)
        %  22 int_rt_rnd_flg:       double, NaN (186)
        %  23 life_cap_rt:          double, NaN (0)
        %  24 loc_amt:              double, NaN (0)
        %  25 max_clm_amt:          double, NaN (0)
        %  26 mnthly_set_aside:     double, NaN (0)
        %  27 mrgn_rt:              double, NaN (0)
        %  28 pmt_md:               double, NaN (0)
        %  29 pmt_strt_dt:          date (elapsed), NaN (197587)
        %  30 prop_addr_zip_cd:     double, NaN (261770)
        %  31 prprty_aprsl_vl:      double, NaN (0)
        %  32 rpr_set_aside:        double, NaN (0)
        %  33 rt_typ:               1, 2, 3, NaN (0)
        %  34 shrd_prem:            double, NaN (0)
        %  35 srvc_fee:             double, NaN (0)
        %  36 srvc_fee_set_aside:   double, NaN (0)
        %  37 stp_late_chg:         double, NaN (0)
        %  38 stp_pnlty_int:        double, NaN (0)
        %  39 tot_loan_bal:         double, NaN (0)
        %  40 trmntn_cd:            31, 32, 33, 34, 35, 36, NaN (327360)
        %  41 trmntn_dt:            date (elapsed), NaN (326623)
        %  42 trmntn_typ:           1, 2, 3, 4, 5, 6, NaN (322937)
        %
        %  After method 'compute' has been invoked:
        %  43 fiscal_year:          double, 0 (1)
        %  44 age:                  double, NaN (0 - after removing)
        %  45 gender:               1, 3, NaN (0 - after removing)        

        % load file
        load hecm_stata

        % transform date information
        actn_dt_y(isinf(actn_dt)) = NaN;
        actn_dt_m(isinf(actn_dt)) = NaN;
        actn_dt_d(isinf(actn_dt)) = NaN;
        actn_dt = datenum(actn_dt_y,actn_dt_m,actn_dt_d);
        actn_dt(isnan(actn_dt))   = 0;

        adjstmnt_dt_y(isinf(adjstmnt_dt)) = NaN;
        adjstmnt_dt_m(isinf(adjstmnt_dt)) = NaN;
        adjstmnt_dt_d(isinf(adjstmnt_dt)) = NaN;
        adjstmnt_dt = datenum(adjstmnt_dt_y,adjstmnt_dt_m,adjstmnt_dt_d);
        adjstmnt_dt(isnan(adjstmnt_dt))   = 0;

        borr_dt_of_brth_y(isinf(borr_dt_of_brth)) = NaN;
        borr_dt_of_brth_m(isinf(borr_dt_of_brth)) = NaN;
        borr_dt_of_brth_d(isinf(borr_dt_of_brth)) = NaN;
        borr_dt_of_brth = datenum(borr_dt_of_brth_y,borr_dt_of_brth_m,borr_dt_of_brth_d);
        borr_dt_of_brth(isnan(borr_dt_of_brth))   = 0;

        borr_dt_of_dth_y(isinf(borr_dt_of_dth)) = NaN;
        borr_dt_of_dth_m(isinf(borr_dt_of_dth)) = NaN;
        borr_dt_of_dth_d(isinf(borr_dt_of_dth)) = NaN;
        borr_dt_of_dth = datenum(borr_dt_of_dth_y,borr_dt_of_dth_m,borr_dt_of_dth_d);
        borr_dt_of_dth(isnan(borr_dt_of_dth))   = 0;

        clsng_dt_y(isinf(clsng_dt)) = NaN;
        clsng_dt_m(isinf(clsng_dt)) = NaN;
        clsng_dt_d(isinf(clsng_dt)) = NaN;
        clsng_dt = datenum(clsng_dt_y,clsng_dt_m,clsng_dt_d);
        clsng_dt(isnan(clsng_dt))   = 0;

        coborr_dt_of_brth_y(isinf(coborr_dt_of_brth)) = NaN;
        coborr_dt_of_brth_m(isinf(coborr_dt_of_brth)) = NaN;
        coborr_dt_of_brth_d(isinf(coborr_dt_of_brth)) = NaN;
        coborr_dt_of_brth = datenum(coborr_dt_of_brth_y,coborr_dt_of_brth_m,coborr_dt_of_brth_d);
        coborr_dt_of_brth(isnan(coborr_dt_of_brth))   = 0;

        coborr_dt_of_dth_y(isinf(coborr_dt_of_dth)) = NaN;
        coborr_dt_of_dth_m(isinf(coborr_dt_of_dth)) = NaN;
        coborr_dt_of_dth_d(isinf(coborr_dt_of_dth)) = NaN;
        coborr_dt_of_dth = datenum(coborr_dt_of_dth_y,coborr_dt_of_dth_m,coborr_dt_of_dth_d);
        coborr_dt_of_dth(isnan(coborr_dt_of_dth))   = 0;

        cs_estblshd_dt_y(isinf(cs_estblshd_dt)) = NaN;
        cs_estblshd_dt_m(isinf(cs_estblshd_dt)) = NaN;
        cs_estblshd_dt_d(isinf(cs_estblshd_dt)) = NaN;
        cs_estblshd_dt = datenum(cs_estblshd_dt_y,cs_estblshd_dt_m,cs_estblshd_dt_d);
        cs_estblshd_dt(isnan(cs_estblshd_dt))   = 0;

        db_dt_y(isinf(db_dt)) = NaN;
        db_dt_m(isinf(db_dt)) = NaN;
        db_dt_d(isinf(db_dt)) = NaN;
        db_dt = datenum(db_dt_y,db_dt_m,db_dt_d);
        db_dt(isnan(db_dt))   = 0;

        endrsmnt_dt_y(isinf(endrsmnt_dt)) = NaN;
        endrsmnt_dt_m(isinf(endrsmnt_dt)) = NaN;
        endrsmnt_dt_d(isinf(endrsmnt_dt)) = NaN;
        endrsmnt_dt = datenum(endrsmnt_dt_y,endrsmnt_dt_m,endrsmnt_dt_d);
        endrsmnt_dt(isnan(endrsmnt_dt))   = 0;

        fndng_dt_y(isinf(fndng_dt)) = NaN;
        fndng_dt_m(isinf(fndng_dt)) = NaN;
        fndng_dt_d(isinf(fndng_dt)) = NaN;
        fndng_dt = datenum(fndng_dt_y,fndng_dt_m,fndng_dt_d);
        fndng_dt(isnan(fndng_dt))   = 0;

        pmt_strt_dt_y(isinf(pmt_strt_dt)) = NaN;
        pmt_strt_dt_m(isinf(pmt_strt_dt)) = NaN;
        pmt_strt_dt_d(isinf(pmt_strt_dt)) = NaN;
        pmt_strt_dt = datenum(pmt_strt_dt_y,pmt_strt_dt_m,pmt_strt_dt_d);
        pmt_strt_dt(isnan(pmt_strt_dt))   = 0;

        trmntn_dt_y(isinf(trmntn_dt)) = NaN;
        trmntn_dt_m(isinf(trmntn_dt)) = NaN;
        trmntn_dt_d(isinf(trmntn_dt)) = NaN;
        trmntn_dt = datenum(trmntn_dt_y,trmntn_dt_m,trmntn_dt_d);
        trmntn_dt(isnan(trmntn_dt))   = 0;

        % replace inf mit NaN
        actn_cd(isinf(actn_cd))                     = NaN;
        borrower_age(isinf(borrower_age))           = NaN;
        int_rt_rnd_flg(isinf(int_rt_rnd_flg))       = NaN;
        prop_addr_zip_cd(isinf(prop_addr_zip_cd))   = NaN;
        rt_typ(isinf(rt_typ))                       = 3;
        trmntn_cd(isinf(trmntn_cd))                 = NaN;
        trmntn_typ(isinf(trmntn_typ))               = NaN;
        coborr_gender(coborr_gender==2)             = NaN;
        borr_gender(borr_gender==2)                 = NaN;


        % create dataset with the Statistics Toolbox
        NumObs    = size(actn_cd,1);
        NameObs   = strcat({'case'},num2str((1:NumObs)','%d'));
        obj.data = dataset({nominal(actn_cd,{'Stop Disbursements','Start Disbursements',...
                                                'Bankruptcy/Litigation','Not Assigned / Label Error (21)',...
                                                'Not Assigned / Label Error (23)','Referred for Foreclosure',...
                                                'Payoff','Liquidated/Held for Sale','Liquidated/Third Party Sale/Condemnation',...
                                                'Assigned to HUD','Not Assigned / Label Error (99)'}),'actn_cd'},...
                            {actn_dt,'actn_dt'},...
                            {adjstmnt_dt,'adjstmnt_dt'},...
                            {borr_dt_of_brth,'borr_dt_of_brth'},...
                            {borr_dt_of_dth,'borr_dt_of_dth'},...
                            {nominal(borr_gender,{'female','male'}),...
                                                  'borr_gender'},...
                            {borrower_age,'borrower_age'},...
                            {clsng_dt,'clsng_dt'},...
                            {coborr_dt_of_brth,'coborr_dt_of_brth'},...
                            {coborr_dt_of_dth,'coborr_dt_of_dth'},...
                            {nominal(coborr_gender,{'female','male'}),...
                                                    'coborr_gender'},...
                            {cs_estblshd_dt,'cs_estblshd_dt'},...
                            {nominal(cs_sts,{'Suspense','Pending Endorsement','Endorsed','Assigned to HUD','Canceled Loan','Terminated','Terminated/Assigned'}),...
                                                    'cs_sts'},...             
                            {db_dt,'db_dt'},...
                            {endrsmnt_dt,'endrsmnt_dt'},...
                            {flr_rt,'flr_rt'},...
                            {fndng_dt,'fndng_dt'},...
                            {init_fee,'init_fee'},...
                            {init_prncpl_lmt,'init_prncpl_lmt'},...
                            {int_rt,'int_rt'},...
                            {int_rt_10yr,'int_rt_10yr'},...
                            {nominal(int_rt_rnd_flg,{'Not Assigned / Label Error (1)','N','Y'}),'int_rt_rnd_flg'},...
                            {life_cap_rt,'life_cap_rt'},...
                            {loc_amt,'loc_amt'},...
                            {max_clm_amt,'max_clm_amt'},...
                            {mnthly_set_aside,'mnthly_set_aside'},...
                            {mrgn_rt,'mrgn_rt'},...
                            {pmt_md,'pmt_md'},...
                            {pmt_strt_dt,'pmt_strt_dt'},...
                            {prop_addr_zip_cd,'prop_addr_zip_cd'},...
                            {prprty_aprsl_vl,'prprty_aprsl_vl'},...
                            {rpr_set_aside,'rpr_set_aside'},...
                            {nominal(rt_typ,{'Annual Adjustable','Monthly Adjustable','Fixed'}),...
                                             'rt_typ'},...                
                            {shrd_prem,'shrd_prem'},...
                            {srvc_fee,'srvc_fee'},...
                            {srvc_fee_set_aside,'srvc_fee_set_aside'},...
                            {stp_late_chg,'stp_late_chg'},...
                            {stp_pnlty_int,'stp_pnlty_int'},...
                            {tot_loan_bal,'tot_loan_bal'},...
                            {nominal(trmntn_cd,{'Death','Borrower is moving to a different residence','Borrower paid off loan but is remaining in the property','Any other reason for termination','Reason for termination unknown','Not Assigned / Label Error (36)'}),...
                            'trmntn_cd'},...                
                            {trmntn_dt,'trmntn_dt'},...
                            {nominal(trmntn_typ,{'Claim','Foreclosure','Not Assigned / Label Error (3)','Payoff','Sale','Voluntary'}),...
                                                'trmntn_typ'},...           
                       'ObsNames',NameObs);

        desc = 'HECM case data provided by E.J.Szymanoski (HUD)';
        obj.data = set(obj.data,'Description',desc);

      end % function CASES
      
      %% OTHERS
      function obj = compute(obj)
        % COMPUTE_CASES  Compute additional variables for HECM database
        %   
        % computes fiscal year of origination, age for calculation of principal
        % limit factor, gender for calculation of principal limit factor; removes
        % cases from dataset where entry of observation is obviously erroneous or
        % data is missing
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Additional Codes: Codierungswert (# of Obs)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % These are the additional variables which are created and appended to
        % database
        %
        %  43 fiscal_year:          double, 0 (1)
        %  44 age:                  double, NaN (0 - after removing)
        %  45 gender:               1, 3, NaN (0 - after removing)

            %% compute fiscal year of origination: 'fiscal_year'
             
                    % load fndng_dt as vectorized date
                    date_start_vec = datevec(obj.data.fndng_dt);

                    % compute fiscal_year
                    fiscal_year = zeros(length(obj.data),1);

                        % use HUD HECM fiscal year definition
                        fiscal_year(date_start_vec(:,2) > 9) = ...
                            date_start_vec((date_start_vec(:,2) > 9),1) + 1;

                        fiscal_year(date_start_vec(:,2) <= 9) =...
                            date_start_vec((date_start_vec(:,2) <= 9),1);   

                    % append to database
                    fiscal_year = dataset({fiscal_year,'fiscal_year'});
                    obj.data = [obj.data fiscal_year]; clear fiscal_year

            %% compute age for plf calculation: 'age'

                    % define availability of birth data           
                    borrower_only = ((obj.data.borr_dt_of_brth~=0)...
                                        &(obj.data.coborr_dt_of_brth==0));

                    coborrower_only = ((obj.data.borr_dt_of_brth==0)...
                                        &(obj.data.coborr_dt_of_brth~=0));

                    both = ((obj.data.borr_dt_of_brth~=0)...
                                &(obj.data.coborr_dt_of_brth~=0));

                    none = ((obj.data.borr_dt_of_brth==0)...
                                &(obj.data.coborr_dt_of_brth==0));

                    % assign relevant date according to case of availability;
                    % POLICY: if both are available, then pick the younger one
                    dt_borrower_only   = borrower_only .* obj.data.borr_dt_of_brth;            
                    dt_coborrower_only = coborrower_only .* obj.data.coborr_dt_of_brth;            
                    dt_both = (both .* (obj.data.borr_dt_of_brth <= obj.data.coborr_dt_of_brth) .* obj.data.borr_dt_of_brth)...
                        + (both .* (obj.data.borr_dt_of_brth > obj.data.coborr_dt_of_brth) .* obj.data.coborr_dt_of_brth);

                    % get relevant date at birth in datevec format
                    dt_birth          = dt_borrower_only + dt_coborrower_only + dt_both;
                    dt_birth_vec      = datevec(dt_birth);
                    dt_birth_vec(:,3) = 1;   % date rounding formular according to HUD Handbook

                    % get date at origination in datevec format
                    dt_origination_vec = datevec(obj.data.fndng_dt);

                    % compute age at origination
                    age_case1 = (~none) .* ((dt_birth_vec(:,2) - dt_origination_vec(:,2)) >= 6) .*...
                        (dt_origination_vec(:,1) - dt_birth_vec(:,1));

                    age_case2 = (~none) .* ((dt_birth_vec(:,2) - dt_origination_vec(:,2)) < 6) .*...
                        (dt_origination_vec(:,1) - dt_birth_vec(:,1) + 1);

                    age       = age_case1 + age_case2;
                    age(none) = NaN;

                    % append to database
                    age = dataset({age,'age'});
                    obj.data = [obj.data age]; clear age

            %% compute gender for plf calculation: 'gender'

                    % preallocate gender content arrays
                    gender = cell(size(obj.data,1),1);

                    % get relevant gender if case is a couple
                    pick_borrower = (both & (obj.data.borr_dt_of_brth <= obj.data.coborr_dt_of_brth));
                    pick_coborrower = (both & (obj.data.borr_dt_of_brth > obj.data.coborr_dt_of_brth));

                    % assign gender for all cases
                    gender(borrower_only) = cellstr(obj.data.borr_gender(borrower_only));
                    gender(coborrower_only) = cellstr(obj.data.coborr_gender(coborrower_only));
                    gender(pick_borrower) = cellstr(obj.data.borr_gender(pick_borrower));
                    gender(pick_coborrower) = cellstr(obj.data.coborr_gender(pick_coborrower));
                    gender(none) = cellstr('<undefined>');

                    % convert gender to nominal array
                    female = strcmp({'female'},gender);
                    male = strcmp({'male'},gender);
                    undefined = strcmp({'<undefined>'},gender);

                    gender_data            = zeros(length(obj.data),1);
                    gender_data(female)    = 1;
                    gender_data(undefined) = NaN;
                    gender_data(male)      = 3;        

                    % append to database
                    gender_nominal = dataset({nominal(gender_data,{'female','male'}),'gender'});
                    obj.data = [obj.data gender_nominal]; clear gender_nominal

             %% remove crazy cases from dataset

                    % define criteria

                    %   1) relevant birth date is inconsistent
                    birthinconsistent  = ((dt_birth_vec(:,1) < 1860) & (dt_birth_vec(:,1) > 0));
                    %   2) age at origination is unrealistic
                    ageinconsistent    = (obj.data.age < 62);
                    %   3) gender could not be matched with birth data
                    noinfo             = isundefined(obj.data.gender);
                    %   4) age at origination could  not be computed
                    noage             = isnan(obj.data.age);
                    %   5) gender at origination could  not be computed
                    nogender          = isundefined(obj.data.gender);

                    % remove
                    obj.data(birthinconsistent,:)    = [];

                    upd_ageinconsistent = ageinconsistent(1:length(obj.data),1);
                    obj.data(upd_ageinconsistent,:)  = [];

                    upd_noinfo = noinfo(1:length(obj.data),1);
                    obj.data(upd_noinfo,:)           = [];

                    upd_noage = noage(1:length(obj.data),1);
                    obj.data(upd_noage,:)           = [];

                    upd_nogender = nogender(1:length(obj.data),1);
                    obj.data(upd_nogender,:)        = [];

                    % show results
                    disp('   ');

                    disp('Cases removed: inconsistent birth data:');
                    disp(int2str(sum(birthinconsistent)));
                    disp('   ');

                    disp('Cases removed: inconsistent age at origination:');
                    disp(int2str(sum(upd_ageinconsistent)));
                    disp('   ');

                    disp('Cases removed: gender could not be matched with missing birth data:');
                    disp(int2str(sum(upd_noinfo)));
                    disp('   ');

                    disp('Cases removed: age at origination could not be computed:');
                    disp(int2str(sum(upd_noage)));
                    disp('   ');

                    disp('Cases removed: relevant gender was not found:');
                    disp(int2str(sum(upd_nogender)));
                    disp('   ');
                    
      end % function COMPUTE
      function varargout = history(obj, observation)
        % HISTORY   How does the time line of a specific contract look like?
        %
        % displays lifetime history of a particular loan from the dataset and
        % returns loan history in a cell array {date, date explanation, date status}

            % read variable names of all dataset variables
            VarName = obj.data.Properties.VarNames';
            DateInd = [2;3;4;5;8;9;10;12;14;15;17;29;41];

            % get variable names of date variables    
            VarDateName = cell(length(DateInd),1);
            for i=1:length(VarDateName),
                VarDateName{i}=VarName{DateInd(i)};
            end

            % get meaning of date variables
            VarDateNameExpl = cell(length(DateInd),1);
                VarDateNameExpl{1} = 'Action Code';
                VarDateNameExpl{2} = 'Interest Rate Adjustment';
                VarDateNameExpl{3} = 'Birth Borrower';
                VarDateNameExpl{4} = 'Death Borrower';
                VarDateNameExpl{5} = 'Closing';
                VarDateNameExpl{6} = 'Birth Coborrower';
                VarDateNameExpl{7} = 'Death Coborrower';
                VarDateNameExpl{8} = 'System Registration';
                VarDateNameExpl{9} = 'Debit Date';
                VarDateNameExpl{10} = 'Endorsement';
                VarDateNameExpl{11} = 'Funds disbursed';
                VarDateNameExpl{12} = 'Monthly Payments Begin';
                VarDateNameExpl{13} = 'Termination';

            VarDateStatus = cell(length(DateInd),1);
                VarDateStatus{1} = char(obj.data.actn_cd(observation));
                VarDateStatus{3} = char(obj.data.borr_gender(observation));
                VarDateStatus{5} = char(obj.data.rt_typ(observation));
                VarDateStatus{6} = char(obj.data.coborr_gender(observation));
                VarDateStatus{13} = char(obj.data.trmntn_typ(observation));

            % get all date values of observation
            dateelapsed = zeros(length(VarDateName),1);
            for i=1:length(VarDateName)
                dateelapsed(i) = obj.data.(VarDateName{i})(observation);
            end

            % sort table, clear missing values
            [dates_sorted, index] = sortrows(dateelapsed);
            index(dates_sorted==0)        = [];
            dates_sorted(dates_sorted==0) = [];

            % set lifetime history for output
            varargout = [cellstr(datestr(dates_sorted)) VarDateNameExpl(index) VarDateStatus(index)];
            disp(varargout);
            disp('    ');
            disp('Current Status:');
            disp(char(obj.data.cs_sts(observation)));

      end % function HISTORY
      function varargout = byyear(obj)
        % BYYEAR   Tabulates number of cases originated by fiscal year
        %
        % Shows and returns a tabulation of HECM contracts originated in a
        % given fiscal year based on the first date of funding 'fndng_dt'

            % load fndng_dt as vectorized date
            date_start_vec = datevec(obj.data.fndng_dt);

            % compute fiscal_year
            fiscal_year = zeros(length(obj.data),1);

                % use HUD HECM fiscal year definition
                fiscal_year(date_start_vec(:,2) > 9) = ...
                    date_start_vec((date_start_vec(:,2) > 9),1) + 1;

                fiscal_year(date_start_vec(:,2) <= 9) =...
                    date_start_vec((date_start_vec(:,2) <= 9),1);

            % sort fiscal years
            fiscal_year_unique = sort(unique(fiscal_year));

            % compute volume originated by fiscal year
            volumeperyear = zeros(length(fiscal_year_unique),1);

                for i=1:length(fiscal_year_unique)
                    volumeperyear(i) = sum(fiscal_year==fiscal_year_unique(i));
                end

            % prepare table output
            printout = {[fiscal_year_unique, volumeperyear]};

            % strip of cases where fndng_dt is NaN
            delete = (printout{1}(:,1)==0);
            ndelete = sum(printout{1}(delete,2)); 
            printout{1}(delete,:) = []; 

            % print how many have been removed from count
            disp('   ');
            disp('# of contracts not being included in table due to missing fndng_dt:');
            disp(int2str(sum(ndelete)));

            % define table
            varargout = printout;

      end % function BYYEAR
      
      %% PLOTTING
      function plot_diff(obj)
          % PLOT_DIFF  plotes the duration in days between closing and
          % funding of contracts (when both dates are known) as boxplot;
          % shows quantiles for duration
          
              % compute difference
              diff = obj.data.fndng_dt((obj.data.fndng_dt ~= 0) & (obj.data.clsng_dt ~= 0))...
                      - obj.data.clsng_dt((obj.data.fndng_dt ~= 0) & (obj.data.clsng_dt ~= 0));
              
              % transform in vector date; get days
              diffvec  = datevec(diff);
              diffdays = 365 .* diffvec(:,1) + 30 .* diffvec(:,2) + diffvec(:,3);              
              
              % quantiles
              start = 100*(0:0.01:0.04);
              main  = 100*(0.05:0.1:0.95);
              final = 100*(0.96:0.01:1);
              percentiles = [start, main, final];
              
              quantiles = prctile(diffdays,percentiles);
              [percentiles; quantiles]'
              
              % plot
              boxplot(diffdays);
      
      end % function PLOT_DIFF

   end % methods  
end % classdef
