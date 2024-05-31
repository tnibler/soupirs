function [DATA,instantaneous_volume]=raw_data_sighs_plot()

% Read and plot the 125 Hz txt file obtained from BlueCherry software


%                                  _______Aur�lien Bringard 04.03.2024

% Filename type "447798_10_1_2020_9_46raw.txt"
% Decimal in dots, column 1 = time, column 2 = flow, column 3 = FO2,  
%column 4 = FCO2
 
% ____________________________ User instructions __________________________

% The function load the data and plot raw data
%__________________________________________________________________________


clear;
day=date;
hour=clock;
parent_function='analyse_soupirs_plot.m';
close all

DIRmat=dir('*.txt'); % list the txt files of the directory

% cd './CPET '
cd './CPET_data_soupirs'

filename=('_ramp-test_raw.log')
% DATA=readmatrix('_ramp-test_raw.log');
DATA=readmatrix(filename);
filename_short=[filename(1:end-8)];


% figure(1);clf
% hold on;
% yyaxis left
% plot(DATA(:,1),DATA(:,4),'g')% FCO2
% plot(DATA(:,1),DATA(:,2),'-r')% flow
% ylabel('flow (l/min) and FCO2 (%)')
% yyaxis right
% plot(DATA(:,1),DATA(:,3),'b');% FO2
% ax=axis;% valeurs min et max des 2 axes: [xmin, xmax, ymin, ymax]
% % trait horizontal noir y=0:
% plot([ax(1),ax(2)],[0,0],'-.k');
% legend('FCO2 (instantaneous)','flow (instantaneous)','FO2 (instantaneous)',...
%     'Location','SouthEast');
% title('raw data')
% % title('select the actual rest period for the fit, using the mouse left button');
% xlabel('time (s)')
% ylabel('FO2 (%)');
% hold off;
% saveas(gcf,[filename_short,'_fig-rawdata.png'])


% Calculation of tidal volume (VT) based on integration of instantaneous flow:
instantaneous_volume=cumtrapz(DATA(:,1),DATA(:,2));
localmax_VT=islocalmax(instantaneous_volume(:),'MinProminence',0.1);
localmin_VT=islocalmin(instantaneous_volume(:),'MinProminence',0.1);
detrend_VT=detrend(instantaneous_volume,3); %detrend VT, to take into account drift
% plot(DATA(:,1),Cumulative_Integral_VE,'b',DATA(:,1),detrend_VT,'g' )


figure(2);clf
hold on;
yyaxis left
plot(DATA(:,1),DATA(:,4),'g')% FCO2
plot(DATA(:,1),DATA(:,2),'-r')% flow
plot(DATA(:,1),instantaneous_volume,'-c')% tidal volume
plot(DATA(localmax_VT,1),instantaneous_volume(localmax_VT),'.k')% local max of VT
plot(DATA(localmin_VT,1),instantaneous_volume(localmin_VT),'.k')% local min of VT
ylabel('flow (l/min), volume (L) and FCO2 (%)')
yyaxis right
plot(DATA(:,1),DATA(:,3),'b');% FO2
ax=axis;% valeurs min et max des 2 axes: [xmin, xmax, ymin, ymax]
% trait horizontal noir y=0:
plot([ax(1),ax(2)],[0,0],'-.k');
legend('FCO2 (instantaneous)','flow (instantaneous)','volume (instantaneous)',...
    'end-inspiration volume','end-expiration volume','FO2 (instantaneous)',...
    'Location','SouthEast');
title('raw data + calculated volume')
xlabel('time (s)')
ylabel('FO2 (%)');
hold off;


time_VTmin=DATA(localmin_VT,1);
time_VTmax=DATA(localmax_VT,1);
size(time_VTmin)
size(time_VTmax)

VTmin=instantaneous_volume(localmin_VT);
VTmax=instantaneous_volume(localmax_VT);
% size(VTmin)
% size(VTmax) 

time_VTmin(1)
time_VTmax(1)

size(time_VTmin,1)-size(time_VTmax,1)

% check the difference between expiratory and inspiratory start points: 
switch (size(time_VTmin,1)-size(time_VTmax,1))
 
% more INSPiratory start points than EXPiratory start points:
    case 1
%     case (size(time_VTmin,1)-size(time_VTmax,1))>=0
        'more INSPiratory start points than EXPiratory start points'
        if time_VTmin(1) < time_VTmax(1) %if first peak is a Vmin   
            'first peak is a Vmin, last peak is a Vmin'
            for   i=1:size((VTmax),1)
                Vinsp(i)=VTmax(i)-VTmin(i);
                Vexp(i)=VTmax(i)-VTmin(i+1);
            end
        elseif time_VTmax(1) < time_VTmin(1) % if first peak is a Vmax           
            for   i=1:size((VTmin),1)
                Vexp(i)=VTmax(i)-VTmin(i);
                Vinsp(i)=VTmax(i)-VTmin(i+1);
            end
        end       
  
     % more EXPiratory start points than INSPiratory start points:   
    case -1
%    case (size(time_VTmin,1)-size(time_VTmax,1))<=0
       'more EXPiratory start points than INSPiratory start points'         
       if time_VTmin(1) < time_VTmax(1) %if first peak is a Vmin   
           for   i=1:size((VTmax),1)
               Vinsp(i)=VTmax(i)-VTmin(i);
               Vexp(i)=VTmax(i)-VTmin(i+1);
           end           
       elseif time_VTmax(1) < time_VTmin(1) % if first peak is a Vmax
           'first peak is a Vmax, last peak is a Vmax'
           for   i=1:size((VTmin),1)
               Vexp(i)=VTmax(i)-VTmin(i);
               Texp(i)=time_VTmin(i)-time_VTmax(i);  
                if i==1
                    Vinsp(i)=[NaN];
                    Tinsp(i)=[NaN];
                    Tcycle(i)=2*Texp(i);                   
                elseif i>1 && i<=size((VTmin),1)
                    Vinsp(i)=VTmax(i)-VTmin(i-1);
                    Tinsp(i)=time_VTmax(i)-time_VTmin(i-1);
                    Tcycle(i)=Tinsp(i)+Texp(i);
                else
                    'problem'
                end
           end
           % delete the last time_VTmax point for Figure 3:
           time_VTmax(size((time_VTmax),1))=[];
       end

       % as much EXPiratory start points than INSPiratory start points:
    case 0 
        'as much EXPiratory start points than INSPiratory start points'
        if time_VTmin(1) < time_VTmax(1) %if first peak is a Vmin
            'first peak is a Vmin, last peak is a Vmax'
            for   i=1:(size((VTmin),1))
%                 i
                Vinsp(i)=VTmax(i)-VTmin(i);
                Tinsp(i)=time_VTmax(i)-time_VTmin(i);
                if i<=(size((VTmax),1)-1)
                    Vexp(i)=VTmax(i)-VTmin(i+1);
                    Texp(i)=time_VTmin(i+1)-time_VTmax(i);
                    Tcycle(i)=Tinsp(i)+Texp(i);
                else
                    'no more VTmin point available'
                    Vexp(i)=[NaN];
                    Texp(i)=[NaN];
                    Tcycle(i)=2*Tinsp(i);
                end                
            end
        elseif time_VTmax(1) < time_VTmin(1) % if first peak is a Vmax
            'first peak is a Vmax, last peak is a Vmin'
            for   i=1:(size((VTmin),1))               
                Vexp(i)=VTmax(i)-VTmin(i);
                Texp(i)=time_VTmin(i)-time_VTmax(i);
                % if i<=(size((VTmin),1)-1)
                if i==1
%                     Vinsp(i)=VTmax(i)-VTmin(i-1);
%                     Tinsp(i)=time_VTmax(i)-time_VTmin(i-1);
%                     Tcycle(i)=Tinsp(i-1)+Texp(i);
                    Vinsp(i)=[NaN];
                    Tinsp(i)=[NaN];
                    Tcycle(i)=2*Texp(i);
                    
%                     Tinsp(i)=time_VTmax(i)-time_VTmin(i-1);
%                     Tcycle(i)=Tinsp(i-1)+Texp(i);
                elseif i>1 && i<=size((VTmin),1)
                    % 'no more VTmax point available'
                    % Vinsp(i)=[NaN];
                    % Tinsp(i)=[NaN];
                    % Tcycle(i)=2*Texp(i);
                    Vinsp(i)=VTmax(i)-VTmin(i-1);
                    Tinsp(i)=time_VTmax(i)-time_VTmin(i-1);
                    Tcycle(i)=Tinsp(i)+Texp(i);
                else
                    'problem'
                end
            end
        end
end

Vinsp=Vinsp';
Tinsp=Tinsp';
Vexp=Vexp';
Texp=Texp';
Tcycle=Tcycle';



% figure(3);
% hold on;
% plot(time_VTmin,Tinsp);
% plot(time_VTmax,Texp);
% legend('Tinsp','Texp');
% hold off;



figure(4);clf % flow volume detection
subplot(2,1,1)
hold on;
% yyaxis left
plot(DATA(:,1),DATA(:,2),'-r')% flow
plot(DATA(:,1),instantaneous_volume,'-c')% tidal volume
ylabel('flow (l/min) and volume (L)')
plot(DATA(localmax_VT,1),instantaneous_volume(localmax_VT),'.k')% local max of VT
plot(DATA(localmin_VT,1),instantaneous_volume(localmin_VT),'.k')% local min of VT
% yyaxis right
ax=axis;% valeurs min et max des 2 axes: [xmin, xmax, ymin, ymax]
plot([ax(1),ax(2)],[0,0],'k');% trait horizontal noir y=0
legend('flow (instantaneous)','volume (instantaneous)','end-inspiration volume',...
    'end-expiration volume','Location','NorthWest');
% title('select the actual rest period for the fit, using the mouse left button');
title('flow volume detection')
xlabel('time (s)')
% ylabel('FO2 (%)');
hold off;

subplot(2,1,2)
hold on;
% plot(time_Vin(1:end-1,1),Vin)
plot(time_VTmin(:,1),Vinsp)
plot(time_VTmax(:,1),Vexp)
% plot(difference)
hold off;
xlabel('time (s)');
ylabel('VT (l)');
legend('Vinsp','Vexp','Location','SouthEast');
saveas(gcf,[filename_short,'_fig-VT-detection.png'])



figure(5);clf % visual detection fo sigh
hold on;
% yyaxis left
plot(DATA(:,1),DATA(:,2),'-r')% flow
plot(DATA(:,1),instantaneous_volume,'-c')% tidal volume
ylabel('flow (l/min) and volume (L)')
plot(DATA(localmax_VT,1),instantaneous_volume(localmax_VT),'.k')% local max of VT
plot(DATA(localmin_VT,1),instantaneous_volume(localmin_VT),'.k')% local min of VT
% yyaxis right
ax=axis;% valeurs min et max des 2 axes: [xmin, xmax, ymin, ymax]

% Add the Vinsp values of each breath:
formatSpec = '%.2f';
txt=[num2str(Vinsp,formatSpec)];
text(DATA(localmax_VT,1),instantaneous_volume(localmax_VT),txt);

plot([ax(1),ax(2)],[0,0],'k');% trait horizontal noir y=0
legend('flow (instantaneous)','volume (instantaneous)','end-inspiration volume',...
    'end-expiration volume','Location','NorthWest');
% title('select the actual rest period for the fit, using the mouse left button');
title('visual detection fo sigh')
xlabel('time (s)')
% ylabel('FO2 (%)');
hold off;



clear localmin_FO2 localmax_FO2 localmin_FCO2 localmax_FCO2 localmin_VT localmax_VT...
    ans ax 

save([filename_short,'_BbB-calculated-from-rawdata_soupirs']);

cd ..