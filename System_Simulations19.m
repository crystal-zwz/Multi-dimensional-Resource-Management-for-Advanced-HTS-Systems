% 二维卫星场景2*3
% - BW-POW: Flexible POW and BW allocation. Fixed user assignment to the beams
% - MAP:    Fixed resources per beam.       Flexible beam user assignment.
% - BW-MAP: Flexible BW allocation.         Flexible beam user assignment.
% - SR
% - SR_BF
% - SR_BW
% - SR_BW_BF

%%
close all
clear all
clc
orig_state = warning;
warning off

%% %  folder for storing
formatOut='yyyy_MM_dd_HH_mm';
text_file=['System_sim_batch' char(datetime('now'),formatOut)]; % Data labelling

% Parent folder for storing the numerical results
parentpath=['System simulations ' char(datetime('now'),formatOut)];
mkdir(parentpath)

%% Selection of resource allocation strategies
en.CR=1;
en.CR_BW=1;
en.MAP=1; % Boolean variable. Enables ( with value equal to 1) the optimization of flexible beam-user mapping( with fixed resources) to cope with the traffic demand.
en.BW_MAP=1; % Boolean variable. Enables( with value equal to 1) the joint optimization of bandwidth and beam-user mapping to cope with the traffic demand.
en.BW_POW=1; % Boolean variable. Enables ( with value equal to 1) the joint optimization of bandwidth and power to cope with the traffic demand.

%%    Two Dimesion scenario
% Monte-Carlo parameters
Nsims_v= 1; % 30h,完整的一轮约30分钟

K=6;% Number beams. IT MUST BE EVEN
% Beam modeling- Bessel
R=50; % km,Beam Radius
beam_cross_roll_off=3;  % dB, roll off value at beam-crossover

%% Place beams in the 2-dimesional scenario
Roll_off_bound=-beam_cross_roll_off;
% b-th波束孔径距n-th用户位置的距离
% 对贝塞尔模型中波束半径为R的波束，在排布局中波束的位置使相邻波束之间的距离为2R
da = 7.0594;
lamda = 299792458/20e9;
theta = atan(50/35786);
u0 = pi*da/lamda*sin(theta);
d=1:0.001:2*R;
u=u0*d/R;

% 从b-th波束到中心波束内n-th用户的信道增益
T = 0.9; % 孔径边缘锥度 aperture edge taper
p = 2;  % 折射率 index of refraction
a = ((p+1)*(1-T)/((p+1)*(1-T)+T))* ...
        (2*besselj(1,u)./ u+ ...
        2^(p+1)*factorial(p)*(T/(1-T))* ...
        besselj((p+1),u)./(u.^(p+1)));
G_dB_ref = 10*log10(abs(a).^2);

% 查找第一个小于-3dB的值，p = 2
ind_roll_off=find(G_dB_ref<=Roll_off_bound,1,'first');
% 可以允许的信噪比下降范围
ind_roll_off_f =find(G_dB_ref<=(8.7-15),1,'first');

% 初始服务半径
R_0=d(ind_roll_off); % 50.015
% 灵活调整后的最大服务半径
R_f=d(ind_roll_off_f); % 71.100

Distance_beam=2*R_0;

% 波束中心坐标
Center_Beams=zeros(K,2);
Center_Beams(1:K/2,1)=0;
Center_Beams(1:K/2,2)=Distance_beam*(0:K/2-1);
Center_Beams(K/2+1:K,1)=Distance_beam*sqrt(3)/2;
Center_Beams(K/2+1:K,2)=Distance_beam*(0:K/2-1)+Distance_beam/2;

%% Traffic distributions definition

% 同质流量
alpha_v{1}=ones(1,K); % Scenario 1
label_alpha{1}='HT'; % Scenario 1 label
mkdir(parentpath,label_alpha{1}) % Create folder to store the data

% 宽热点
alpha_v{2}=[1 1 1 5 5 1]; % Scenario 3
label_alpha{2}='WHS';  % Scenario 3 label
mkdir(parentpath,label_alpha{2})  % Create folder to store the data

Nscen=2;

%% System parameters

M=4; % Number of carrier per colour ( 4-colour scheme in the 2-dimensional scenario)
Delta_W=  1/(2*M); % Portion of the carrier bandwidth

total_band=500e6; % Total Bandwidth
% Traffic requested per user in Mbps,
% normalized to the avalaible bandwidth : 25 Mbps/ 500 MHz
Req_user_ref= 25e6/(total_band);
% log2(1+ SNR_eff_beam) ,log(1+10^1.35)/log(2) with SNR_eff_beam the average SNR for
% a uniform distribution of user within a beam ( approx 13.5 dB )
SE_av=4.5271;
% Average spectral effiency with a two-Colour scheme 6.8e3/500/6
SE_beam_uniform=SE_av/2;
max_user_per_beam = round(SE_beam_uniform/Req_user_ref);

% Number of user to achieve the capacity per beam
Nuser_beam=  SE_beam_uniform./Req_user_ref;

% h_total小时中随时间变化的用户数di
% 24h
% Nuser_tot = [140,95,60,40,25,17,45,80,130,165,182,205, ...
%     225,245,276,285,290,299,295,285,265,240,200,180];
Nuser_tot = [10,80,150,220 ...
             290,360,430,500];

%% Satellite Parameters

% Maximum saturation power
Max_P_sat= 200*K/6;
% Max Saturation downlink power per Sat in Watts ( beam)
Max_Pb_Sat = 2*Max_P_sat/K;
% Reference value of Saturation downlink power per Sat in Watts
Ref_Pb_Sat = Max_P_sat/K; % 200 W for 6 beams
% Max Power per HPA
Max_Pamp=Ref_Pb_Sat*2*2;

% Satellite repeater output loss (dB)
Tx_sat_repeater_output_loss_dB = 2;
% Satellite repeater antenna loss (dB)
Tx_sat_repeater_antenna_loss_dB = 0.05;
% Downlink Polarization Loss (dB
DL_polarization_loss_dB = 0.2;
% Maximum satellite antenna gain (dBi)
G_max_Sat= 52;

% Terminal Antenna Diameter (m)
terminal_antenna_diameter = 0.6;
% Terminal Antenna Efficiency (<1)
terminal_antenna_efficiency = 0.65;
% Losses due to terminal depointing (dB)
terminal_depointing_losses = 0.5;
% Terminal component losses (dB)
terminal_component_loss_dB = 0;
% Terminal component Noise Temperature (K)
terminal_component_T = 310;
% Noise Figure of the terminal LNB (dB)
terminal_LNA_NF_dB = 2;


symbol_rate= total_band;

% Average free space losse (dB)
FSL= 210;
% Average atmospheric losses (dB)
L_atm=0.4296;
% Ground Noise Temperature (K)
ground_Tground = 45;
% Average clear sky noise temperature at the receiver(K)
T_sky_av=28.4082;
% Average cloud noise temperature at the receiver(K)
T_cloud_av=0.6712;

%%% Constants
% Terminal temperature (K)
terminal_default_T = 290;
% Speed of light (m/s)
speed_of_light = 299792458;
% Boltzmann Constant (J/K)
boltzmann_constant = 1.3806503e-23;

%%% Input parameters used for the generation of the provided files, do not change
% Downlink Frequency Hz
DL_freq = 20e9;

%% Link Budget

% Compute the antenna temperature
T_ta_rx =T_sky_av+T_cloud_av;
% Compute the total noise temperature for each point on the grid
T_tot_rx    = T_ta_rx+ground_Tground+(10^(terminal_component_loss_dB/10)-1)*terminal_component_T+(10^(terminal_LNA_NF_dB/10)-1)*terminal_default_T/(10^(-terminal_component_loss_dB/10));

% Compute the maximum antenna gain for the user
terminal_antenna_Gmax = 10*log10(terminal_antenna_efficiency.*(pi*terminal_antenna_diameter*DL_freq/speed_of_light).^2);


% Compute the G/kT (in dB) for each grid point
GkT  = terminal_antenna_Gmax- FSL - DL_polarization_loss_dB - terminal_depointing_losses-10*log10(T_tot_rx)-10*log10(boltzmann_constant);
% Compute a refernce C/N factor (in dB) without considering the antenna pattern (assumic isotropic antenna) and without the transmit power
L_CN0_dB = -Tx_sat_repeater_output_loss_dB - Tx_sat_repeater_antenna_loss_dB + GkT - L_atm;

% Compute a refernce SNR factor taking into account the symbol rate
symbol_rate_dB = 10*log10(symbol_rate);
L_SNR_dB = L_CN0_dB - symbol_rate_dB;

%% Solver variables
Max_time_opt_second=10; % s, Maxium simulation time for the second-step process. One value for each simulated scenario
optionsGA.PoP= 2000*2; % Population for the genetic algorithm
optionsGA.Elite_num=ceil(0.1*optionsGA.PoP);% Number of elite members for the genetic algorithm
optionsGA.pmut=0.1; % Mutation probability for the genetic algorithm
optionsGA.Max_Gen=5000; % Maximum number of generations for the genetic algorithm


%% Batch simulations
for ind_scen=Nscen  % Select Scenario
    Nsims = Nsims_v;
    % Variables that stores the simulated data
    Nuser_beam_c=cell(1,Nsims);
    Req_user_c=cell(1,Nsims);
    R_off_CR_c=cell(1,Nsims);
    R_off_CR_PAA_c=cell(1,Nsims);
    R_off_MAP_c=cell(1,Nsims);
    R_off_BW_CR_c=cell(1,Nsims);
    R_off_BW_CR_PAA_c=cell(1,Nsims);
    R_off_BW_MAP_c=cell(1,Nsims);
    R_off_BW_POW_c=cell(1,Nsims);
    Assig_MAP_c=cell(1,Nsims);
    Assig_CR_c=cell(1,Nsims);
    x_ga_paa_m = zeros(Nsims,2*K);
    x_ga_paa_m1 = zeros(Nsims,2*K);
    Assig_BW_CR_c=cell(1,Nsims);
    Assig_BW_MAP_c=cell(1,Nsims);
    users_locations_c=cell(1,Nsims);
    snr_car_c=cell(1,Nsims);
    user_beams_c=cell(1,Nsims);
    Bandwidth_Allo_BW_CR=cell(1,Nsims);
    Bandwidth_Allo_BW_MAP=cell(1,Nsims);
    Power_Allo_BW_POW=cell(1,Nsims);
    Bandwidth_Allo_BW_POW=cell(1,Nsims);

    res_IterNum_Av = 0;
    res_NQU_Av=0;
    res_NU_Av=0;
    res_OffRate_Av=0;
    res_MinUserRate_Av=0;

    sprintf('Start of Scenario %i',ind_scen) % Display to keep track of the simulations
    tic
    for indsims=1:Nsims % Monte-Carlo simulations for the selected scenario and batch size.

        sprintf('Monte-Carlo simulation %i of %i',indsims,Nsims) % Display to keep track of the simulations

        %%   Number of user per beam
        aux_alpha_v= alpha_v{ind_scen}; % Values to model the traffic. Dirichlet distribution

        % Obtain random number following a Dirichlet distribution
        n=1; % Number of generated scenarios
        % 从 gamma 分布中生成一个随机数数组
        r = gamrnd(repmat(aux_alpha_v,n,1),1,n,K); % Generate n vectors of longitude K following a gamma distribution Gamma(alpha_i,1)
        Dir_rand = r ./ repmat(sum(r,2),1,K); % Normalization to obtain random numbers following a Dirichlet distribution.

        % Obtain a integer number of users per beam
        Nuser_beam=round( max(Nuser_tot)*Dir_rand);
        N=sum(Nuser_beam);
        % 避免round使数量不一致
        Nuser_tot(Nuser_tot == max(Nuser_tot)) = N;
        % Auxiliar variable 辅助变量
        cum_Nu_sim_index=cumsum(Nuser_beam);

        %%  User location generation   %
        % User are randomly placed within the beam radius
        % 生成一个用户最多时刻的随机分布
        user_beams=zeros(1,N);   % Vector that indicates the domminant beam for each user
        users_b_index=cell(K,1); % Cell with the user indexes for each beam
        users_b_loc=cell(K,1);
        for i=1:K
            x0_beam=Center_Beams(i,1);
            y0_beam=Center_Beams(i,2);

            if Nuser_beam(i)~=0
                % Generate random user locations within a beam
                t = 2*pi*rand(Nuser_beam(i),1);
                r = R*sqrt(rand(Nuser_beam(i),1)); % Beam with radius R from Bessel modeling
                % Obtain user position
                x = x0_beam + r.*cos(t);
                y = y0_beam + r.*sin(t);
                % Generate auxiliar variables
                switch i
                    case 1
                        users_b_index{i}=1:cum_Nu_sim_index(i);
                    otherwise
                        users_b_index{i}=cum_Nu_sim_index(i-1)+1:cum_Nu_sim_index(i);
                end
                user_beams(users_b_index{i})=i;
                users_b_loc{i}=[x y];
            end
        end
        users_locations=cell2mat(users_b_loc);

        %%  Obtain channel values

        %%% Compute distances between user locations and beam centers for Bessel modeling
        distance_All=zeros(N,K);
        for i=1:K
            distance_All(:,i)= vecnorm(users_locations-Center_Beams(i,:),2,2);
        end

        %%% Aplay Bessel modeling
        u=u0*distance_All/R;
        indZ=find(~distance_All);

        g = ((p+1)*(1-T)/((p+1)*(1-T)+T))* ...
            (2*besselj(1,u)./ u+ ...
            2^(p+1)*factorial(p)*(T/(1-T))* ...
            besselj((p+1),u)./(u.^(p+1)));
        G_dB = 10*log10(abs(g).^2);
        
        G_dB(indZ)=zeros(1,length(indZ));

        %%% Magnitude Channel matrix, normalized respect to the Noise with bandwidth W
        Gamma_dB_wo_P=L_SNR_dB+(G_dB+G_max_Sat); %  in dB
        Gamma_wo_P=10.^(Gamma_dB_wo_P/10); % in natural

        % Carrier SNR for uniform resource allocation
        snr_car= Gamma_wo_P*Ref_Pb_Sat/0.5; % Carrier SNR
        snr_car_unfilter=snr_car;
        SNR_car=10*log10(snr_car); % 272*6
        Th_car= 8.7000; % SNR threshold to obtain a C/I=24 dB or higher.
        ind_snr = false(N,K);
        for rol = 1:N
            ind_snr(rol,:)=SNR_car(rol,:)<Th_car;
            if min(ind_snr(rol,:))==1
                [~,idmax] = max(SNR_car(rol,:));
                ind_snr(rol,idmax) = 0;
            end
        end
        snr_car(ind_snr)=0;   % Filter SNR values
        
        %% Structure with CR
        scen_data.CenterBeams_CR = Center_Beams;
        scen_data.CenterBeams_BWCR = Center_Beams;
        
        if en.CR
            f1 = figure('name',strcat('SR: Monte-Carlo simulation ',num2str(indsims)," of ",num2str(Nsims)));
            t1 = tiledlayout(2,4);
            % 'flow'
        end

        if en.CR_BW
            f2 = figure('name',strcat('BW-SR: Monte-Carlo simulation ',num2str(indsims)," of ",num2str(Nsims)));
            t2 = tiledlayout(2,4);
        end

        if en.BW_MAP
            f3 = figure('name',strcat('BW_MAP: Monte-Carlo simulation ',num2str(indsims)," of ",num2str(Nsims)));
            t3 = tiledlayout(2,4);
        end

        h_total = length(Nuser_tot);
        %% 一天h_total小时的流量变化
        for h = 1:h_total
        % for h = 18 % test
            %%  Auxiliar varibles for the resource managment  %
            N_h = Nuser_tot(h);
            % 根据每小时用户数，从最大用户数的分布中随机选取对应数量的点
            N_id = sort(randperm(max(Nuser_tot),N_h)); 
            cum_Nu_sim_index_h1 = [0,cum_Nu_sim_index];
            for k = 1:K
                users_b_index_h1{k} =N_id(((cum_Nu_sim_index_h1(k)+1)<=N_id)&(N_id<=cum_Nu_sim_index_h1(k+1)));
                Nuser_beam_h(k) =length(users_b_index_h1{k});
            end
            cum_Nu_sim_index_h2=cumsum(Nuser_beam_h);
            cum_Nu_sim_index_h3 = [0,cum_Nu_sim_index_h2];
            for k = 1:K
                users_b_index_h{k} = (cum_Nu_sim_index_h3(k)+1):cum_Nu_sim_index_h3(k+1);
            end
             
            user_beams_h = user_beams(N_id);
            users_locations_h = users_locations(N_id,:);
            distance_All_h = distance_All(N_id,:);
            Gamma_wo_P_h = Gamma_wo_P(N_id,:);
            snr_car_unfilter_h = snr_car_unfilter(N_id,:);
            snr_car_h = snr_car(N_id,:);
            
            ind_sims_h = h_total*(indsims-1)+h;

            % Requested traffic per user
            Req_user=Req_user_ref*ones(N,1);
            Req_user_h=Req_user(N_id);
            Req_user_c{ind_sims_h}=Req_user_h+0;

            Total_Req=sum(Req_user_h); % Total requested traffic

            % Traffic requested per beam
            Req_b=zeros(1,K);

            for i=1:K
                users_index=users_b_index{i};
                if ~isempty(users_index)
                    Req_b(i)=sum(Req_user(users_index));
                end
            end

            %% Structure with the simulated scenario
            scen_data.K=K;
            scen_data.M=M;
            scen_data.N=N_h;
            scen_data.Delta_W=Delta_W;
            scen_data.Nuser_beam=Nuser_beam_h;
            scen_data.user_beams=user_beams_h;
            scen_data.users_b_index=users_b_index_h;
            scen_data.Max_P_sat=Max_P_sat;
            scen_data.Ref_Pb_Sat=Ref_Pb_Sat;
            scen_data.Max_Pamp=Max_Pamp;
            scen_data.aux_Cn=(0.5/Ref_Pb_Sat)*snr_car_unfilter_h;
            scen_data.Req_beam=Req_b;
            scen_data.Req_user=Req_user_h;
            scen_data.Gamma_wo_P=Gamma_wo_P_h;
            scen_data.snr_car= snr_car_h;
            scen_data.Max_time_opt_second=Max_time_opt_second;
            scen_data.SE_av=SE_av;

            scen_data.R_0=R_0;
            scen_data.Center_Beams=Center_Beams;
            scen_data.R_f=R_f;
            scen_data.distance_All=distance_All_h;
            scen_data.max_user_per_beam=max_user_per_beam;
            scen_data.users_locations=users_locations_h;
            scen_data.L_SNR_dB=L_SNR_dB;

            % Create a "ResourceAssignment" object  with the simulated data
            Sim_object= ResourceAssignment19(scen_data);

            %%              Resource assignment                %%
            %% SR
            if en.CR
                sprintf('Solution: Flexible RangeCover')
                [R_off_CR,R_off_CR_PAA,x_ga_paa,Assig_CR,iterationNum_CR,CenterBeams_CR,Assig_UM_FlexB,Radius_Beams_CR] = FlexibleCoverRange(Sim_object);
                scen_data.CenterBeams_CR=CenterBeams_CR;

                % Store user rates
                R_off_CR_c{ind_sims_h}=R_off_CR+0;
                R_off_CR_PAA_c{ind_sims_h}=R_off_CR_PAA+0;
                % Store beam-user mapping
                Assig_CR_c{ind_sims_h}=Assig_CR';
                % Store beam gain and power assignment
                x_ga_paa_m(ind_sims_h,:) = x_ga_paa;

                % Obtain the Normalized Quadratic unment demand (NQU)
                NQU_CR=sum(  ( Req_user_h- R_off_CR).^2  )*(total_band^2)/( N_h*(total_band*Req_user_ref)^2 );
                NQU_CR_PAA=sum(  ( Req_user_h- R_off_CR_PAA).^2  )*(total_band^2)/( N_h*(total_band*Req_user_ref)^2 );               
                % Obtain the Normalized Unment demand (NU)
                NU_CR=sum(  ( Req_user_h- R_off_CR)  )/(SE_beam_uniform*K);
                NU_CR_PAA=sum(  ( Req_user_h- R_off_CR_PAA)  )/(SE_beam_uniform*K);                
                % Obtain the total offered rate, Mbps
                Offer_rate_CR=sum( R_off_CR)*total_band/1e6;
                Offer_rate_CR_PAA=sum( R_off_CR_PAA)*total_band/1e6;
                % Obtain the minimum user rate offered rate, Mbps
                Min_User_rate_CR=min( R_off_CR)*total_band/1e6;
                Min_User_rate_CR_PAA=min( R_off_CR_PAA)*total_band/1e6;
                
                % 波束服务范围时变
                figure(f1);
                nexttile
                for ind_beam=1:K                   
                    if Assig_CR(ind_beam)> 0
                        %% 绘制用户位置和波束服务范围
                        users_index=find(Assig_UM_FlexB==ind_beam);
                        x_CR = users_locations_h(users_index,1);
                        y_CR = users_locations_h(users_index,2);
                        color = ["#574BB4","#C0321A","#629C35","#DD7C4F","#911eb4","#6F6F6F"];

                        center_CR = CenterBeams_CR(ind_beam,:);
                        radius_CR = Radius_Beams_CR(ind_beam);
                        pos = [center_CR(1)-radius_CR, ...
                            center_CR(2)-radius_CR, ...
                            2*radius_CR, 2*radius_CR]; 
                        r1 = rectangle('Position',pos,'Curvature',[1 1],'LineStyle','-');
                        r1.EdgeColor = color(ind_beam);
                        r1.LineWidth = 1.2;
                        hold on

                        scatter(x_CR,y_CR,...
                            20,'x','MarkerEdgeColor',color(ind_beam),'LineWidth',1.5);
                        hold on
                        title(strcat(num2str(5+h-1),":00 ",'No.',num2str(iterationNum_CR)))
                        ax = gca;
                        ax.TitleHorizontalAlignment = 'left';

                        ylim tight
                        xlim tight
                        axis equal

                        %% 绘制波束方向图
                        df=-R:1:R;
                        uf=2.07123*df/R;
                        theta = linspace(0,2*pi,length(df));
                        [Theta,U] = meshgrid(theta,uf);
                        gp = ((x_ga_paa(ind_beam)+1)*(1-T)/((x_ga_paa(ind_beam)+1)*(1-T)+T))* ...
                            (2*besselj(1,U)./ U+ ...
                            2^(x_ga_paa(ind_beam)+1)*factorial(x_ga_paa(ind_beam))*(T/(1-T))* ...
                            besselj((x_ga_paa(ind_beam)+1),U)./(U.^(x_ga_paa(ind_beam)+1)));
                        G_dB_f = 10*log10(abs(gp).^2);
                        
                        Gamma_wo_P_CR_f= L_SNR_dB+G_dB_f+G_max_Sat+10*log10(x_ga_paa(ind_beam+6));
                        s = surf(U/2.0712*50.*cos(Theta),U/2.0712*50.*sin(Theta),Gamma_wo_P_CR_f);
                        % 使 z 轴方向的一个数据单位的长度等于 x 轴方向和 y 轴方向的n个数据单位的长度
                        daspect([50/2.0712 50/2.0712 1])
                        % FaceAlpha_p = [0.33,0.25,0.1];
                        % s.FaceAlpha = FaceAlpha_p(x_ga_paa(ind_beam));
                        s.FaceAlpha = 0.2;
                        s.EdgeColor = 'none';
                        sx = get(s,'xdata');    
                        sy = get(s,'ydata');
                        % sz = get(s,'zdata');
                        set(s,'xdata',sx+center_CR(1));    
                        set(s,'ydata',sy+center_CR(2));
                        % set(s,'zdata',sz+x_ga_paa(ind_beam+6));
                        hold on
                    end
                    % view(75,35)
                    view(2)

                end

            else
                NQU_CR=NaN;
                NU_CR=NaN;
                Offer_rate_CR=NaN;
                Min_User_rate_CR=NaN;

                NQU_CR_PAA=NaN;
                NU_CR_PAA=NaN;
                Offer_rate_CR_PAA=NaN;
                Min_User_rate_CR_PAA=NaN;

                iterationNum_CR = NaN;
            end

            %% Flexible Beam-user mapping with fixed resources. Two-step optimization process.
            if en.MAP
                sprintf('Solution: Flexible beam-user mapping with fixed resources')
                [R_off_MAP,Assig_MAP]=FixResFlexibleMapping(Sim_object);
                % Store user rates
                R_off_MAP_c{ind_sims_h}=R_off_MAP+0;
                % Store beam-user mapping
                Assig_MAP_c{ind_sims_h}=Assig_MAP;

                % Obtain the Normalized Quadratic unment demand (NQU)
                NQU_MAP=sum(  ( Req_user_h- R_off_MAP).^2  )*(total_band^2)/( N_h*(total_band*Req_user_ref)^2 );
                % Obtain the Normalized Unment demand (NU)
                NU_MAP=sum(  ( Req_user_h- R_off_MAP)  )/(SE_beam_uniform*K);
                % Obtain the total offered rate, Mbps
                Offer_rate_MAP=sum( R_off_MAP)*total_band/1e6;
                % Obtain the minimum user rate offered rate, Mbps
                Min_User_rate_MAP=min( R_off_MAP)*total_band/1e6;
            else
                NQU_MAP=NaN;
                NU_MAP=NaN;
                Offer_rate_MAP=NaN;
                Min_User_rate_MAP=NaN;
            end

            %% SR_BW
            if en.CR_BW
                sprintf('Solution: Flexible bandwidth and RangeCover')
                [R_off_BW_CR,R_off_BW_CR_PAA,x_ga_paa1,Assig_BW_CR,M_BW_CR,iterationNum_CRBW,CenterBeams_BWCR,Assig_UM_FixedRes,Radius_Beams_BWCR] = FlexBandwidthFlexibleCoverRange(Sim_object);
                scen_data.CenterBeams_BWCR=CenterBeams_BWCR;
                % Store user rates
                R_off_BW_CR_c{ind_sims_h}=R_off_BW_CR+0;
                R_off_BW_CR_PAA_c{ind_sims_h}=R_off_BW_CR_PAA+0;
                % Store beam-user mapping
                Assig_BW_CR_c{ind_sims_h}=Assig_BW_CR';
                % Store bandwidth allocation
                Bandwidth_Allo_BW_CR{ind_sims_h}=M_BW_CR';
                % Store beam gain and power assignment
                x_ga_paa_m1(ind_sims_h,:) = x_ga_paa1;

                % Obtain the Normalized Quadratic unment demand (NQU)
                NQU_BW_CR=sum(  ( Req_user_h- R_off_BW_CR).^2  )*(total_band^2)/( N_h*(total_band*Req_user_ref)^2 );
                NQU_BW_CR_PAA=sum(  ( Req_user_h- R_off_BW_CR_PAA).^2  )*(total_band^2)/( N_h*(total_band*Req_user_ref)^2 );
                % Obtain the Normalized Unment demand (NU)
                NU_BW_CR=sum(  ( Req_user_h- R_off_BW_CR)  )/(SE_beam_uniform*K);
                NU_BW_CR_PAA=sum(  ( Req_user_h- R_off_BW_CR_PAA)  )/(SE_beam_uniform*K);
                % Obtain the total offered rate, Mbps
                Offer_rate_BW_CR=sum( R_off_BW_CR)*total_band/1e6;
                Offer_rate_BW_CR_PAA=sum( R_off_BW_CR_PAA)*total_band/1e6;
                % Obtain the minimum user rate offered rate, Mbps
                Min_User_rate_BW_CR=min( R_off_BW_CR)*total_band/1e6;
                Min_User_rate_BW_CR_PAA=min( R_off_BW_CR_PAA)*total_band/1e6;
                
                figure(f2);
                nexttile
                for ind_beam=1:K
                    if Assig_BW_CR(ind_beam)> 0
                        users_index=find(Assig_UM_FixedRes==ind_beam);
                        x_CR = users_locations_h(users_index,1);
                        y_CR = users_locations_h(users_index,2);
                        color = ["#574BB4","#C0321A","#629C35","#DD7C4F","#911eb4","#6F6F6F"];

                        center_CR = CenterBeams_BWCR(ind_beam,:);
                        radius_CR = Radius_Beams_BWCR(ind_beam);
                        pos = [center_CR(1)-radius_CR, ...
                            center_CR(2)-radius_CR, ...
                            2*radius_CR, 2*radius_CR];
                        r1 = rectangle('Position',pos,'Curvature',[1 1],'LineStyle','-');
                        r1.EdgeColor = color(ind_beam);
                        r1.LineWidth = 1.2;
                        hold on

                        scatter(x_CR,y_CR,...
                            20,'x','MarkerEdgeColor',color(ind_beam),'LineWidth',1.5);
                        hold on
                        title(strcat(num2str(5+h-1),":00 ",'No.',num2str(iterationNum_CRBW)))
                        ax = gca;
                        ax.TitleHorizontalAlignment = 'left';

                        ylim tight
                        xlim tight
                        axis equal

                        %% 绘制波束方向图
                        df=-R:1:R;
                        uf=2.07123*df/R;
                        theta = linspace(0,2*pi,length(df));
                        [Theta,U] = meshgrid(theta,uf);
                        gp = ((x_ga_paa1(ind_beam)+1)*(1-T)/((x_ga_paa1(ind_beam)+1)*(1-T)+T))* ...
                            (2*besselj(1,U)./ U+ ...
                            2^(x_ga_paa1(ind_beam)+1)*factorial(x_ga_paa1(ind_beam))*(T/(1-T))* ...
                            besselj((x_ga_paa1(ind_beam)+1),U)./(U.^(x_ga_paa1(ind_beam)+1)));
                        G_dB_f = 10*log10(abs(gp).^2);
                        Gamma_wo_P_CR_f= L_SNR_dB+G_dB_f+G_max_Sat+10*log10(x_ga_paa1(ind_beam+6));
                        s = surf(U/2.0712*50.*cos(Theta),U/2.0712*50.*sin(Theta),Gamma_wo_P_CR_f);
                        % 使 z 轴方向的一个数据单位的长度等于 x 轴方向和 y 轴方向的n个数据单位的长度
                        daspect([50/2.0712 50/2.0712 1])
                        % FaceAlpha_p = [0.33,0.25,0.1];
                        % s.FaceAlpha = FaceAlpha_p(x_ga_paa1(ind_beam));
                        s.FaceAlpha = 0.2;
                        s.EdgeColor = 'none';
                        sx = get(s,'xdata');    
                        sy = get(s,'ydata');
                        % sz = get(s,'zdata');
                        set(s,'xdata',sx+center_CR(1));    
                        set(s,'ydata',sy+center_CR(2));
                        % set(s,'zdata',sz+x_ga_paa1(ind_beam+6));
                        hold on
                    end
                    % view(75,35)
                    view(2)
                end

            else
                NQU_BW_CR=NaN;
                NU_BW_CR=NaN;
                Offer_rate_BW_CR=NaN;
                Min_User_rate_BW_CR=NaN;
                iterationNum_CRBW = NaN;

                NQU_BW_CR_PAA=NaN;
                NU_BW_CR_PAA=NaN;
                Offer_rate_BW_CR_PAA=NaN;
                Min_User_rate_BW_CR_PAA=NaN;
            end

            %% Flexible Bandwidth and Beam-user mapping. Two-step optimization process.
            if en.BW_MAP
                sprintf('Solution: Flexible bandwidth and beam-user mapping')
                [R_off_BW_MAP,Assig_BW_MAP,M_BW_MAP,Assig_UM_FixedRes_MAP]=FlexBandwidthFlexMapping(Sim_object);
                % Store user rates
                R_off_BW_MAP_c{ind_sims_h}=R_off_BW_MAP+0;
                % Store beam-user mapping
                Assig_BW_MAP_c{ind_sims_h}=Assig_BW_MAP;
                % Store bandwidth allocation
                Bandwidth_Allo_BW_MAP{ind_sims_h}=M_BW_MAP';

                % Obtain the Normalized Quadratic unment demand (NQU)
                NQU_BW_MAP=sum(  ( Req_user_h- R_off_BW_MAP).^2  )*(total_band^2)/( N_h*(total_band*Req_user_ref)^2 );
                % Obtain the Normalized Unment demand (NU)
                NU_BW_MAP=sum(  ( Req_user_h- R_off_BW_MAP)  )/(SE_beam_uniform*K);
                % Obtain the total offered rate, Mbps
                Offer_rate_BW_MAP=sum( R_off_BW_MAP)*total_band/1e6;
                % Obtain the minimum user rate offered rate, Mbps
                Min_User_rate_BW_MAP=min( R_off_BW_MAP)*total_band/1e6;

                 figure(f3);
                nexttile
                for ind_beam=1:K
                    if Assig_BW_MAP(ind_beam)> 0
                        users_index=find(Assig_UM_FixedRes_MAP==ind_beam);
                        x_CR = users_locations_h(users_index,1);
                        y_CR = users_locations_h(users_index,2);
                        color = ["#574BB4","#C0321A","#629C35","#DD7C4F","#911eb4","#6F6F6F"];

                        center_CR = Center_Beams(ind_beam,:);
                        radius_CR = R_0;
                        pos = [center_CR(1)-radius_CR, ...
                            center_CR(2)-radius_CR, ...
                            2*radius_CR, 2*radius_CR];
                        r1 = rectangle('Position',pos,'Curvature',[1 1],'LineStyle','-');
                        r1.EdgeColor = color(ind_beam);
                        r1.LineWidth = 1.2;
                        hold on

                        scatter(x_CR,y_CR,...
                            20,'x','MarkerEdgeColor',color(ind_beam),'LineWidth',1.5);
                        hold on
                        title(strcat(num2str(5+h-1),":00 "))
                        ax = gca;
                        ax.TitleHorizontalAlignment = 'left';

                        ylim tight
                        xlim tight
                        axis equal
                    end
                end
            else
                NQU_BW_MAP=NaN;
                NU_BW_MAP=NaN;
                Offer_rate_BW_MAP=NaN;
                Min_User_rate_BW_MAP=NaN;
            end

            % Flexible Bandwidth and Power. Genetic Algorithm
            if en.BW_POW
                sprintf('Solution: Flexible bandwidth and power')

                % Default values if other strategies with fixed mapping are not simualted.
                P_POW=Ref_Pb_Sat*ones(1,K);
                M_BW=M*ones(K,1);

                [R_off_BW_POW,M_BW_POW,P_BW_POW]=FlexBandwidthPower(Sim_object,P_POW,M_BW,optionsGA);
                % Store user rates
                R_off_BW_POW_c{ind_sims_h} =R_off_BW_POW+0;
                % Store bandwidth allocation
                Bandwidth_Allo_BW_POW{ind_sims_h}=M_BW_POW;
                % Store power allocation
                Power_Allo_BW_POW{ind_sims_h}=P_BW_POW;

                % Obtain the Normalized Quadratic unment demand (NQU)
                NQU_BW_POW=sum(  ( Req_user_h- R_off_BW_POW).^2  )*(total_band^2)/( N_h*(total_band*Req_user_ref)^2 );
                % Obtain the Normalized Unment demand (NU)
                NU_BW_POW=sum(  ( Req_user_h- R_off_BW_POW)  )/(SE_beam_uniform*K);
                % Obtain the total offered rate, Mbps
                Offer_rate_BW_POW=sum( R_off_BW_POW)*total_band/1e6;
                % Obtain the minimum user rate offered rate, Mbps
                Min_User_rate_BW_POW=min( R_off_BW_POW)*total_band/1e6;
            else
                NQU_BW_POW=NaN;
                NU_BW_POW=NaN;
                Offer_rate_BW_POW=NaN;
                Min_User_rate_BW_POW=NaN;
            end

            % Store basic data to replicate the scenario
            user_beams_c{ind_sims_h}=user_beams_h;
            Nuser_beam_c{ind_sims_h}=Nuser_beam_h;
            users_locations_c{ind_sims_h}=users_locations_h;

            % Display Normalized Quadratic unment demand (NQU)

            sprintf('Simulation Results:')
            res_IterNum = [ iterationNum_CR iterationNum_CRBW];
            res_NQU = [ NQU_BW_POW NQU_MAP NQU_BW_MAP NQU_CR NQU_CR_PAA NQU_BW_CR NQU_BW_CR_PAA]';
            res_NU =  [ NU_BW_POW NU_MAP NU_BW_MAP NU_CR NU_CR_PAA NU_BW_CR NU_BW_CR_PAA]';
            res_OffRate = [ Offer_rate_BW_POW Offer_rate_MAP Offer_rate_BW_MAP ...
                Offer_rate_CR Offer_rate_CR_PAA Offer_rate_BW_CR Offer_rate_BW_CR_PAA]';
            res_MinUserRate = [ Min_User_rate_BW_POW Min_User_rate_MAP Min_User_rate_BW_MAP ...
                Min_User_rate_CR Min_User_rate_CR_PAA Min_User_rate_BW_CR Min_User_rate_BW_CR_PAA]';

            Table = table(res_NQU,res_NU,res_OffRate,res_MinUserRate,'VariableNames', ...
                {'NQU','NU','Offered Rate [Mps]','Min. User Rate [Mbps]'}, ...
                'RowName',{'BW-POW','MAP','BW-MAP','SR','SR-BF','SR-BW','SR-BW-BF'});
            disp(Table)

            res_IterNum_Av = res_IterNum_Av + res_IterNum;
            res_NQU_Av=res_NQU_Av+res_NQU;
            res_NU_Av=res_NU_Av+res_NU;
            res_OffRate_Av=res_OffRate_Av+res_OffRate;
            res_MinUserRate_Av=res_MinUserRate_Av+res_MinUserRate;
            toc

            sprintf("End of %i o'clock in Scenario %i",h+4 ,ind_scen) % Display to keep track of the simulations
        
        end % loop hours

        if en.CR
            figure(f1);
            cb = colorbar;
            cb.Layout.Tile = 'east';
            cb.Label.String = 'beam gain [dB]';
            xlabel(t1,'x [km]')
            ylabel(t1,'y [km]')
            t1.TileSpacing = 'tight';
            t1.Padding = 'tight';
        end

        if en.CR_BW
            figure(f2);
            cb = colorbar;
            cb.Layout.Tile = 'east';
            cb.Label.String = 'beam gain [dB]';
            xlabel(t2,'x [km]')
            ylabel(t2,'y [km]')
            t2.TileSpacing = 'tight';
            t2.Padding = 'tight';
        end

        if en.BW_MAP
            figure(f3);
            xlabel(t3,'x [km]')
            ylabel(t3,'y [km]')
            t3.TileSpacing = 'tight';
            t3.Padding = 'tight';
        end
        
    end % loop Nsims

    sprintf('Batch simulations Results:')
    T = table(res_NQU_Av/Nsims/h,res_NU_Av/Nsims/h,res_OffRate_Av/Nsims/h,res_MinUserRate_Av/Nsims/h, ...
        'VariableNames',{'NQU','NU','Offered Rate [Mps]','Min. User Rate [Mbps]'}, ...
        'RowName',{'BW-POW','MAP','BW-MAP','SR','SR-BF','SR-BW','SR-BW-BF'});
    disp(T)
    IterNum_Av = res_IterNum_Av/Nsims/h;
    disp(IterNum_Av) % 200M,1.0394


    sprintf('Storing Batch simulation results')
    save([ parentpath  '\' label_alpha{ind_scen} '\'  text_file '_alpha_' label_alpha{ind_scen} '.mat'],...
        'Nsims_v','M','K','Delta_W','Nuser_beam_c','Req_user_c','R_off_CR_c','R_off_CR_PAA_c','R_off_MAP_c', ...
        'R_off_BW_CR_c','R_off_BW_CR_PAA_c','R_off_BW_MAP_c','Bandwidth_Allo_BW_CR','Bandwidth_Allo_BW_MAP', ...
        'Power_Allo_BW_POW','Bandwidth_Allo_BW_POW','user_beams_c','users_locations_c','Assig_BW_CR_c', ...
        'Assig_MAP_c','Assig_CR_c','Assig_BW_MAP_c','R_off_BW_POW_c','IterNum_Av','h_total','x_ga_paa_m1','x_ga_paa_m')


end % loop Nscen


warning(orig_state)
