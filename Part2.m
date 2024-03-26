clear
%% Radical SR1 specs  
WB=2.217;           %wheelbase (m)
radius=50;          %skidpan radius (m)
g=9.81;             %gravity
Cl=0.27;             %lift coefficient
m=588;             %car mass (kg)
rtol_wd=0.48;       %right to left weight distribution
ltor_wd=1-rtol_wd;  %left to right weight distribution
CG=0.338;           %cg height (m)
Rch=0.237;          %roll center height (m)
RearTW=1.343;       %Rear track width (m)
FrontTW=1.345;      %Front track width (m)
A=1.04;             %frontal area (m^2);
ro=1.225;           %air density (kg/m^3)
AB=0.45;            %Aero balance
v_vector=0:0.1:33;  %velocity vector
Kfroll=60000;       %front roll rate Nm/deg
Krroll=60000;       %rear roll rate Nm/deg
WD=0;       %rear to front weight distribution
j=0;
for v=0:0.1:33 % increasing speed by 0.1 m/s at every loop
    LatG=(v^2/(radius*g));       %Lateral g calculations velocity ^2 / (skid pan radius * g) (N)
    j=j+1;
    LatG_v(j)=LatG;             %Lateral G vector
        %% Downforce calculations
    Fl=0.5*Cl*A*ro*v^2;             % downforce formula (N)
    Fl_front=Fl*AB;                     % downforce on front axis (N)
    Fl_rear=Fl-Fl_front;                   %downforce on rear axis (N)

    i=0;
for WD=0.45:0.01:0.50 
    i=i+1;
    b=WB*(1-WD);        %distance of front axle from CG (m)   
    a=WB*WD;            %distance of rear axle from CG (m)   

    %% Geometric Weight transfer calculations
    FGWT(j,i)=(((m*(a/WB))*(LatG)*Rch)/(FrontTW))*g;     %Front axis Lateral weight transfer calclulation -> ( mass * (a/wheelbase) * LatG * roll center height) / Front track width (N)
    RGWT(j,i)=(((m*(b/WB))*(LatG)*Rch)/(RearTW))*g;     %Rear axis Lateral weight transfer calculation -> (mass * (b/wheelbase) * LatG * roll center height) / Rear track width (N)
    
    %% Elastic weight tranfer calculations
    FEWT(j,i)=((((m-40)*(LatG)*(CG-Rch))*g/(FrontTW))*(Kfroll/(Kfroll+Krroll))); %elastic weight transfer on front axis (N)
    REWT(j,i)=((((m-40)*(LatG)*(CG-Rch))*g/(RearTW))*(Krroll/(Krroll+Kfroll))); %elastic weight transfer on rear axis (N)
    
    
    %% Total weight transfer caclulations
    FWT(j,i)=FGWT(j,i)+FEWT(j,i);                  %Front Weight Transfer (N)
    RWT(j,i)=RGWT(j,i)+REWT(j,i);                  %Rear Weight Transfer (N)
    
    %% Dynamic Load at each tyre
    FzFL(j,i)=((((m*WD*rtol_wd)*g)-FWT(j,i))+(Fl_front))/1000;            %Front left tyre dynamic load  (kN)
    FzFR(j,i)=((((m*WD*ltor_wd)*g)+FWT(j,i))+(Fl_front))/1000;             %Front right tyre dynamic load (kN)
    FzRL(j,i)=((((m*WD*rtol_wd)*g)-RWT(j,i))+(Fl_rear))/1000;            %Rear left tyre dynamic load  (kN)
    FzRR(j,i)=((((m*WD*ltor_wd)*g)+RWT(j,i))+(Fl_rear))/1000;            %Rear right tyre dynamic load (kN)

    %% Magic Formula

    % % % %Front Tyre Coefficient Parameters % % %
    ay1=10; ay2=1300; ay3=1100; ay4=10; ay6=0; ay7=-2; 
    Cy=1.4;

    %Lateral Force Calculations
    
    % % % Front Right Tyre Lateral Force Calculation % % %
    Dy=((FzFR(j,i))*(ay1*FzFR(j,i)))+ay2; 
    BCDy = ay3*sind(atan(FzFR(j,i)/ay4)*2); 
    By=BCDy/(Cy*Dy); 
    Ey=(ay6*(FzFR(j,i)))+ay7;
    alphaf=rad2deg((FzFR(j,i)*v^2)/(BCDy*g*radius)); %front slip angle  
    FyFR(j,i)=Dy*sind(Cy*atand(By*(alphaf)-Ey*(By*(alphaf)-atand(By*(alphaf)))));
    
    %Front Left Tyre Lateral Force Calculation % % %
    Dy=FzFL(j,i)*(ay1.*FzFL(j,i)+ay2); 
    BCDy = ay3*sind(atan(FzFL(j,i)/ay4)*2); 
    By=BCDy/(Cy*Dy); 
    Ey=(ay6*FzFL(j,i))+ay7; 
    alphaf=rad2deg((FzFL(j,i)*v^2)/(BCDy*g*radius)); %front slip angle  
    FyFL(j,i)=Dy*sind(Cy*atand(By*(alphaf)-Ey*(By*(alphaf)-atand(By*(alphaf)))));
    
    % % % Rear tyre coeffcient parameters % % %
    ay1=0; ay2=1100; ay3=1100; ay4=10; ay6=0; ay7=-2; 
    Cy=1.4;

    % % % Rear Right Tyre Lateral Force Calculation % % %
    Dy=FzRR(j,i)*(ay1*FzRR(j,i))+ay2; 
    BCDy = ay3*sind(atan((FzRR(j,i)))/ay4)*2; 
    By=BCDy/(Cy*Dy); 
    Ey=(ay6*(FzRR(j,i)))+ay7; 
    alphar=rad2deg((FzRR(j,i)*v^2)/(BCDy*g*radius)); %rear slip angle 
    FyRR(j,i)=Dy*sind(Cy*atand(By*(alphar)-Ey*(By*(alphar)-atand(By*(alphar)))));
    
    % % % Rear Left Tyre Lateral Force Calculation % % %
    Dy=FzRL(j,i)*(ay1*FzRL(j,i)+ay2); 
    BCDy = ay3*sind(atan(FzRL(j,i)/ay4)*2); 
    By=BCDy/(Cy*Dy); 
    Ey=(ay6*FzRL(j,i))+ay7;  
    alphar=rad2deg((FzRL(j,i)*v^2)/(BCDy*g*radius)); %rear slip angle 
    FyRL(j,i)=Dy*sind(Cy*atand(By*(alphar)-Ey*(By*(alphar)-atand(By*(alphar)))));
    
    %% Yaw Moment Calculation   
    YawMom(j,i)=(a*(FyFL(j,i)+FyFR(j,i)))-(b*(FyRL(j,i)+FyRR(j,i)));
    lgnd{i} = ['WD = ', num2str(WD*100),'%'];
end
end
%% Plot Dynamic Load on each tyre
figure('Position', [0 50 1000 900])
t = tiledlayout(2,2);
title(t,'Dynamic Load on each corner')
xlabel(t,'Lateral Acceleration (g)')
ylabel (t,'Vertical Load Fz (kN)')
nexttile
plot(LatG_v,FzFL)
title('Front Left')
legend(lgnd,'location','best') 
grid on
nexttile
plot(LatG_v,FzFR)
title('Front Right')
legend(lgnd,'location','best') 
grid on
nexttile
plot(LatG_v,FzRL)
title('Rear Left')
legend(lgnd,'location','best') 
grid on
nexttile
plot(LatG_v,FzRR)
legend(lgnd,'location','best') 
grid on
title('Rear Right')

%% Yaw Moment Plot
figure('Position', [1100 50 700 700])
plot(LatG_v,YawMom, 'LineWidth',1.5)
title('Yaw Moment Plot')
xlabel('Lateral Acceleration (g)')
ylabel ('Yaw Moment (Nm')
line(xlim(), [0,0], 'LineWidth', 1, 'Color', 'k');
grid on;
legend(lgnd,'location','best') 
%% Plot generated Lateral Force on each tyre
figure('Position', [0 50 1000 900])
t = tiledlayout(2,2);
title(t,'Lateral Force on each corner')
xlabel(t,'Lateral Acceleration (g)')
ylabel (t,'Lateral Force (N)')
nexttile
plot(LatG_v,FyFL)
title('Front Left')
legend(lgnd,'location','best') 
grid on
nexttile
plot(LatG_v,FyFR)
title('Front Right')
legend(lgnd,'location','best') 
grid on
nexttile
plot(LatG_v,FyRL)
title('Rear Left')
legend(lgnd,'location','best') 
grid on
nexttile
plot(LatG_v,FyRR)
legend(lgnd,'location','best') 
grid on
title('Rear Right')

