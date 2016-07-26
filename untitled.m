function varargout = untitled(varargin)
% UNTITLED MATLAB code for untitled.fig
%      UNTITLED, by itself, creates a new UNTITLED or raises the existing
%      singleton*.
%
%      H = UNTITLED returns the handle to a new UNTITLED or the handle to
%      the existing singleton*.
%
%      UNTITLED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNTITLED.M with the given input arguments.
%
%      UNTITLED('Property','Value',...) creates a new UNTITLED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help untitled

% Last Modified by GUIDE v2.5 07-Nov-2013 13:19:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @untitled_OpeningFcn, ...
                   'gui_OutputFcn',  @untitled_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before untitled is made visible.
function untitled_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to untitled (see VARARGIN)

% Choose default command line output for untitled
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using untitled.
%if strcmp(get(hObject,'Visible'),'off')
%    plot(rand(5));
%end

% UIWAIT makes untitled wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = untitled_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Simulation based on the model in Gross, Pugh, Burns 2012
global R0
global k_R
global k_R2
global k_Rfraction
global k_E
global k_E2
global k_Efraction
global nu_RE
global D_Ca
global D_cG
global beta_dark
global n_cyc
global K_cyc
global alpha_max
global cG_dark
global Ca_dark
global J_cGdark_f2
global J_cGdark
global FA2B
global J_satex
global K_ex
global delta_beta4
global alpha_dark
global I_dark
global beta_idv
global n_cG %n_cG=3 in Gross, Pugh, Burns 2012
n_cG=3;
%global K_cG
%K_cG=20; %microM, TODO from Kolesnikov et al 2011

%Set parameters
popup_sel_index = get(handles.popupmenu2, 'Value');
switch popup_sel_index
    case 1
        Ca_feedback=1
    case 2
        Ca_feedback=0
end

R0=str2double(get(handles.edit_Rstar,'String')) %Number of Rhodopsin isomerizations

%ROS geometry, see Parameter values section
delta_interdiscal=0.031 %micron, interdiscal distance
L_OS=22 %micron (micrometer), rod outer segment length
A_OS=1.7 %micron^2, rod outer segment cross sectional area 

N_av=6.022141 * 10^(23); %Avogadro's number
e_q=1.6021766 * 10^(-19); %magnitude of electron charge, in Coulombs
F=N_av*e_q %Faraday constant

%See Table 2

%{
k_R= 25 %1/sec, Rate of R* activation
k_E= 5 %1/sec, WT, Rate of G*-E* deactivation
%k_E= 12.5; %RGS9-ox
nu_RE= 300; %1/sec, Max rate of G*-E* activation per R*
%}
k_R=str2double(get(handles.edit_k_R,'String'))
k_Rfraction=str2double(get(handles.edit_kR1fraction,'String'))
k_R2=str2double(get(handles.edit_kR2,'String'))
k_E=str2double(get(handles.edit_k_E,'String'))
k_Efraction=str2double(get(handles.edit_kE1fraction,'String'))
k_E2=str2double(get(handles.edit_k_E2,'String'))
nu_RE=str2double(get(handles.edit_nu_RE,'String'))

%{
beta_idv= 43; %1/sec, Rate of cGMP hydrolysis per G*-E*
beta_dark= 4.1; %1/sec, Rate of spontaneous cGMP hydrolysis
%}
beta_idv=str2double(get(handles.edit_beta_idv,'String'))
beta_dark=str2double(get(handles.edit_beta_dark,'String'))

D_Ca= 2 %micron^2/s, see Parameter values section
D_cG= 40 %micron^2/s, Longitudinal diffusion coefficient of cGMP

%{
f_Ca= 0.12; %Fraction of current carried by calcium
B_Ca= 50; %Calcium buffer capacity
%}
f_Ca=str2double(get(handles.edit_f_Ca,'String'))
B_Ca=str2double(get(handles.edit_B_Ca,'String'))

%{
alpha_dark= 16.7; %microM/sec, Dark cGMP synthesis rate (Where is this used in the model??? In Eq. (9))
n_cyc= 1.5; %Hill coefficient for Ca2+ dependence of cGMP synthesis
K_cyc= 0.080; %microM, K1/2 for Ca2+ dependence of cGMP synthesis
alpha_max= 150; %microM/sec, Max rate of cGMP synthesis
%}
alpha_dark=str2double(get(handles.edit_alpha_dark,'String'))
alpha_max=str2double(get(handles.edit_alpha_max,'String'))
n_cyc=str2double(get(handles.edit_n_cyc,'String'))
K_cyc=str2double(get(handles.edit_K_cyc,'String'))

%{
cG_dark=4.1 %microM
Ca_dark=0.320 %microM
K_ex= 1.1; %microM, K1/2 for NCKX activation
J_satex= 0.21; %pA/micron, Maximum NCKX current
%Table 1
I_dark=16.7; %pA, GCAPS-/-
%I_dark=15.4; %pA, GCAPS-/-RGS9-ox
%Gross, Pugh and Burns did not provide the value for J_cGdark:(
%so we need to deduce or infer its value
%Hypothesis 1: J_cGdark may be calculated from Eq. (8) by setting r(t)=0
%J_cGdark=I_dark/L_OS - J_satex*Ca_dark/(Ca_dark+K_ex) %0.7118 pA/micron
%Hypothesis 2: J_cGdark may be calculated from Eq. (2) by setting square bracketed
%term to 0.
J_cGdark=2/f_Ca*J_satex*Ca_dark/(Ca_dark+K_ex) %0.7887 pA/micron
%cG_dark=alpha_max/(1+(Ca_dark/K_cyc)^n_cyc)/beta_dark %Gives 4.0650, close
%to 4.1
%cG_dark=alpha_dark/beta_dark %Gives 4.0732, close to 4.1
%}
cG_dark=str2double(get(handles.edit_cG_dark,'String'))
Ca_dark=str2double(get(handles.edit_Ca_dark,'String'))
%TEST changing cooperativity
%cG_dark=alpha_max/(1+(Ca_dark/K_cyc)^n_cyc)/beta_dark %Gives 4.0650, close
J_cGdark=str2double(get(handles.edit_J_cGdark,'String'))
K_ex=str2double(get(handles.edit_K_ex,'String'))
J_satex=str2double(get(handles.edit_J_satex,'String'))
I_dark=str2double(get(handles.edit_I_dark,'String'))
%TODO Calculate like Dr. Burns did
%J_cGdark=I_dark/(L_OS*(1+f_Ca/2));
%J_satex=f_Ca/2*J_cGdark*(Ca_dark+K_ex)/Ca_dark;

%J_cGdark=2/f_Ca*J_satex*Ca_dark/(Ca_dark+K_ex) %0.7887 pA/micron
%I_dark=L_OS*(J_cGdark+J_satex*Ca_dark/(Ca_dark+K_ex)) %18.3932 pA

%Collect factors (See Eq. (2)), and put conversion factor from mol/micron^3 to microMolar (microM or micromols/liter)
%which is 10^(21), to try to reduce multiplications during pde solving.
%FA2B=10^(-12)/(F*A_OS*0.5*B_Ca)*10^(21); %First factor converts pA to A.
FA2B=10^9/(F*A_OS*0.5*B_Ca)
%
J_cGdark_f2=J_cGdark*f_Ca*0.5;
%See Eq. (6)
delta_beta4=delta_interdiscal*beta_idv/4

numtimedivs=str2double(get(handles.edit_timediv,'String'))
%Plot the time course of E*
%tau_E=1/k_E;
TimeAfterFlash=str2double(get(handles.edit_timeafter,'String'))
%dt=tau_E*10/numtimedivs
dt=TimeAfterFlash/numtimedivs
T=0:dt:TimeAfterFlash; %Time points to evaluate solutions
axes(handles.axes1);
cla;
plot(T,E(T),'g','Linewidth',5);
title('Time course of active PDE');
xlabel('Time t (s)');
ylabel('E*(t)');
%hold on
%hr=plot(T,R0*exp(-k_R*T),'Linewidth',5);
%set(hr,'Color',[1 0.5 0.2]);
%hold off

%Plot the spatiotemporal course of cG and Ca for x>x0
%x0=L_OS/2;
x0=L_OS*get(handles.slider1,'Value')
RightHalf=L_OS-x0;
numrightdiv=str2double(get(handles.edit_x0rightdiv,'String'))
dx_right=RightHalf/numrightdiv
Xright=x0:dx_right:L_OS;
[u1_right,u2_right]=pde_call_right(Xright,T,Ca_feedback);
axes(handles.axes3);
cla;
%figure
rotate3d on; %TODO How to turn rotate3d off for other panels (e.g. 2D plots)?
surf(Xright,T,u1_right)
title('cGMP, cG(x,t)')
xlabel('Distance x (microns)')
ylabel('Time t (s)')
zlabel('cGMP (microMolar)');

axes(handles.axes4);
cla;
%figure
surf(Xright,T,u2_right)
title('Calcium, Ca(x,t)')
xlabel('Distance x (microns)')
ylabel('Time t (s)')
zlabel('Ca (microMolar)');

%Plot the spatiotemporal course of cG and Ca for x<x0
%x0=L_OS/2;
LeftHalf=x0;
numleftdiv=str2double(get(handles.edit_x0leftdiv,'String'))
dx_left=x0/numleftdiv
Xleft=0:dx_left:x0;
[u1_left,u2_left]=pde_call_left(Xleft,T,Ca_feedback);
axes(handles.axes3);
hold on
%cla;
%figure
surf(Xleft,T,u1_left)
%title('cG(x,t)')
%xlabel('Distance x')
%ylabel('Time t')
hold off

axes(handles.axes4);
hold on
%cla;
%figure
surf(Xleft,T,u2_left)
%title('Ca(x,t)')
%xlabel('Distance x')
%ylabel('Time t')
hold off
%rotate3d off;

%Test discontinuity
axes(handles.axes5);
%hcGMP=figure;
cGt=str2num(get(handles.edit_timepoint_cG,'String'))
plot([u1_left(cGt,:),u1_right(cGt,:)]/cG_dark,[Xleft,Xright]); %The first index in the array
%refers to time.
title('[cGMP]/[cGMP]D');
%Plot closed-form solution (Eq. 12)
%figure(hcGMP);
axes(handles.axes5);
hold on
plot(Rogue_SteadyState_cG([Xleft,Xright],x0),[Xleft,Xright],'--r');
hold off
%dlmwrite('Rogue_sim_pdepe.dat',transpose([[Xleft,Xright];[u1_left(cGt,:),u1_right(cGt,:)]/cG_dark]),'delimiter',' ');

axes(handles.axes2);
cla;
rt=LightResponse_PointSink(u1_left,u2_left,u1_right,u2_right,dx_left,dx_right);
%rt=alpha_minus_beta_cG(u1_left,u2_left,u1_right,u2_right,dx_left,dx_right); %test, d[cG]/dt, spatially integrated (spatial average?)
%{
SS=size(T);
Z=zeros(SS)
for i=1:SS(2)
    Z(i,:)=sin(2*3.14159*i/10)/L_OS*ones(size(T));
    %Z(i,:)=i/L_OS*ones(size(T));
end
rt=LightResponse_PointSink(Z,Z,Z,Z,dx_left,dx_right);
%}
%dlmwrite('test.dat',transpose([T;E(T);rt]),'delimiter',' ');
%hburns=figure %Temporary
%plot(T(1:length(T)-1),rt,'g','LineWidth',5);
plot(T,rt,'g','LineWidth',5);
title('Light response (Difference between Dark current and current during the response)');
xlabel('Time t (s)')
ylabel('r(t) (pA)')
%dlmwrite('pdepe_GCAPsmmRGS9ox.dat',transpose([T;rt]),'delimiter',' ');
%dlmwrite('Feedback_tau15ms_ncyc3.dat',transpose([T;rt]),'delimiter',' ');
%dlmwrite('Feedback_response.dat',transpose([T(1:length(T)-1);rt]),'delimiter',' ');
%dlmwrite('NoFeedback_response.dat',transpose([T;rt]),'delimiter',' ');


%hold on
%hkraft=figure
%KraftData=load('kraft.dat');
%plot(KraftData(:,1),-2*[KraftData(:,2),KraftData(:,3),KraftData(:,4),KraftData(:,5),KraftData(:,6),KraftData(:,7)])

%figure(hburns);
%hold on
%BurnsData=load('SPR_burns.dat');
%plot(BurnsData(:,1)-1.025,[BurnsData(:,2),BurnsData(:,3),BurnsData(:,4),BurnsData(:,5)])
%plot(BurnsData(:,1)-1.025,[BurnsData(:,2),BurnsData(:,3)])
%plot(BurnsData(:,1)-1.025,[BurnsData(:,4),BurnsData(:,5)])
%hold off

%Record cGMP and Calcium time course
%dlmwrite('test.dat',transpose([T ; transpose(u1_right(:,1)) ; transpose(u2_right(:,1))]),'delimiter',' ');

%JcG=arrayfun(@JcG_func,u1_right(:,1));
%Jex=arrayfun(@Jex_func,u2_right(:,1));
%dlmwrite('test.dat',transpose([T ; transpose(JcG) ; transpose(Jex)]),'delimiter',' ');

%Special solution, light strikes everywhere along the outer segment
check_alldisks = get(handles.checkbox_alldisks, 'Value');
if check_alldisks
    numalldiv=numleftdiv+numrightdiv
    dx=L_OS/numalldiv
    Xall=0:dx:L_OS;
    [u1_all,u2_all]=pde_call_alldisks(Xall,T,Ca_feedback);

    axes(handles.axes2);
    %figure(hkraft);
    rt_all=LightResponse_alldisks(u1_all,u2_all,dx);
    %dlmwrite('test.dat',transpose([T;E(T);rt_all]),'delimiter',' ');
    hold on
    plot(T,rt_all,'g','LineWidth',5);
    hold off

    %hold on
    %KraftData=load('kraft.dat');
    %plot(KraftData(:,1),-2*[KraftData(:,2),KraftData(:,3),KraftData(:,4),KraftData(:,5),KraftData(:,6),KraftData(:,7)])
    %hold off

    axes(handles.axes5);
    hold on
    plot(u1_all(cGt,:)/cG_dark,Xall,'g');
    hold off
%{
    %TODO Check that differential equation is being satisfied
    figure
    hold on
    xi=length(u1_all(1,:))-10;
    expkRt=exp(-k_R*T(2:end));
    %disp test, u1_all(2:end,xi)
    expkR2t=exp(-k_R2*T(2:end));
    expkEt=exp(-k_E*T(2:end));
    expkE2t=exp(-k_E2*T(2:end));
    E1=k_Rfraction*R0*nu_RE/(k_E-k_R)*(expkRt-expkEt)+(1-k_Rfraction)*R0*nu_RE/(k_E-k_R2)*(expkR2t-expkEt);
    %disp test, E1
    E2=k_Rfraction*R0*nu_RE/(k_E2-k_R)*(expkRt-expkE2t)+(1-k_Rfraction)*R0*nu_RE/(k_E2-k_R2)*(expkR2t-expkE2t);
    %Transpose!
    RHS=alpha_max./(1+(u2_all(2:end,xi)/K_cyc).^n_cyc)-(beta_dark + 0.5*beta_idv*(k_Efraction*E1+(1-k_Efraction)*E2) )'.*u1_all(2:end,xi);
    plot(T(2:end),diff(u1_all(:,xi))/dt,'--r');
    plot(T(2:end),RHS); %Differential equation indeed satisfied!
    %alpha_max./(1+(u2_all(2:end,xi)/K_cyc).^n_cyc)-b(beta_dark + 0.5*beta_idv*E1 ).*u1_all(2:end,xi));%+D_cG*(u1_left(2:end,xi+1)-2*u1_left(2:end,xi)+u1_left(2:end,xi-1))/(dx_left*dx_left);
    hold off
%}
end
%{
%Test if Neumann BC at x0 is being satisfied
figure
hold on
plot(T,delta_beta4*E(T).*(u1_left(:,end))');
%plot(T,delta_beta4*E(T).*(u1_right(:,1))');
plot(T,-D_cG*(u1_left(:,end)-u1_left(:,end-1))/dx_left,'--r'); %Looks good
%plot(T,D_cG*(u1_right(:,2)-u1_right(:,1))/dx_right,'--g');
hold off

%Test if the zero flux BC at 0 and L are being satisfied (I also checked
%these for the alldisks case.
figure
hold on
plot(T,-D_cG*(u1_left(:,2)-u1_left(:,1))/dx_left,'--r');
plot(T,D_cG*(u1_right(:,end)-u1_right(:,end-1))/dx_right,'--g');
hold off

%Wow! Differential equations also satisfied in the SPR model!
%For cGMP
figure
hold on
xi=length(u1_left(1,:))-5;
plot(T(2:end),diff(u1_left(:,xi))/dt,'--r');
%plot(T(2:end),(u1_left(2:end,xi)),'--g');
%plot(T(2:end),alpha_dark-beta_dark*u1_left(2:end,xi) + D_cG*(u1_left(2:end,xi+1)-2*u1_left(2:end,xi)+u1_left(2:end,xi-1))/(dx_left*dx_left));
plot(T(2:end),alpha_max./(1+(u2_left(2:end,xi)/K_cyc).^n_cyc)-beta_dark*u1_left(2:end,xi)+D_cG*(u1_left(2:end,xi+1)-2*u1_left(2:end,xi)+u1_left(2:end,xi-1))/(dx_left*dx_left));
%plot(T(2:end),1./(1+((u1_left(2:end,xi))/K_cyc).^n_cyc));
hold off

%For Ca
figure
hold on
xi=length(u1_left(1,:))-5;
plot(T(2:end),diff(u2_left(:,xi))/dt,'--r');
plot(T(2:end),FA2B*(J_cGdark_f2*(u1_left(2:end,xi)/cG_dark).^n_cG-J_satex*u2_left(2:end,xi)./(u2_left(2:end,xi)+K_ex))+D_Ca*(u2_left(2:end,xi+1)-2*u2_left(2:end,xi)+u2_left(2:end,xi-1))/(dx_left*dx_left));
hold off
%}

%Equation (12)
%Solution to the steady-state profile of cG(x), with R*=1
%Uses k_E. Does not use k_E2!!!
function y=Rogue_SteadyState_cG(x,x0);
global k_E
global nu_RE
%global cG_dark
global D_cG
global beta_dark
global delta_beta4
%WARNING! Does not use k_E2
y=1-exp(-abs(x-x0)*(beta_dark/D_cG)^0.5)/(1+4*k_E*(beta_dark*D_cG)^0.5/(nu_RE*delta_beta4*4));

%Active PDE, E*, Eq. (7b)
function y=E(t);
global R0
global k_R
global k_R2
global k_Rfraction
global k_E
global k_E2
global k_Efraction
global nu_RE

expkRt=exp(-k_R*t);
expkR2t=exp(-k_R2*t);
expkEt=exp(-k_E*t);
expkE2t=exp(-k_E2*t);
if 1
    E1=k_Rfraction*R0*nu_RE/(k_E-k_R)*(expkRt-expkEt)+(1-k_Rfraction)*R0*nu_RE/(k_E-k_R2)*(expkR2t-expkEt);
else
    %Special combination
    %The intermediate R2 is created at rate k_R, while R is lost
    %at rate k_R.
    %The following expression should be symmetric in k_R and k_R2
    E1=k_R*(R0*nu_RE/(k_E-k_R)*(expkRt-expkEt)-R0*nu_RE/(k_E-k_R2)*(expkR2t-expkEt))/(k_R2-k_R) + R0*nu_RE/(k_E-k_R)*(expkRt-expkEt);
    %To make sure
    %k_Efraction=1;
end
E2=k_Rfraction*R0*nu_RE/(k_E2-k_R)*(expkRt-expkE2t)+(1-k_Rfraction)*R0*nu_RE/(k_E2-k_R2)*(expkR2t-expkE2t);
y=k_Efraction*E1+(1-k_Efraction)*E2;
%y=k_Rfraction*R0*nu_RE/(k_E-k_R)*(exp(-k_R*t)-exp(-k_E*t))+(1-k_Rfraction)*R0*nu_RE/(k_E-k_R2)*(exp(-k_R2*t)-exp(-k_E*t));

%Light response, r(t), Eq. (8)
%(I have checked that the numerical integration is robust to the number of
%subdivisions.)
function r=LightResponse_PointSink(cG_left,Ca_left,cG_right,Ca_right,dx_left,dx_right);
global I_dark
%Do the interval [0,x0]
%Calculate currents for each spatiotemporal point
JcG=arrayfun(@JcG_func,cG_left);
Jex=arrayfun(@Jex_func,Ca_left);
%Integrate over the length of the outer segment for each t
sum_left=trapz(transpose(JcG+Jex));
%Do the interval [x0,L]
%Calculate currents for each spatiotemporal point
JcG=arrayfun(@JcG_func,cG_right);
Jex=arrayfun(@Jex_func,Ca_right);
%Integrate over the length of the outer segment for each t
sum_right=trapz(transpose(JcG+Jex));
r=I_dark-(sum_left*dx_left+sum_right*dx_right);
%r=I_dark-2*(sum_left*dx_left);

function r=LightResponse_alldisks(cG,Ca,dx);
global I_dark
%Do the interval [0,L]
%Calculate currents for each spatiotemporal point
JcG=arrayfun(@JcG_func,cG);
Jex=arrayfun(@Jex_func,Ca);
%Integrate over the length of the outer segment for each t
sum=trapz(transpose(JcG+Jex));
r=I_dark-(sum*dx);

%Eq. 4
function y=JcG_func(x);
global J_cGdark
global cG_dark
global n_cG
y=J_cGdark*(x/cG_dark)^n_cG; %See Eq. (4)
%y=0;

%Eq. 5
function y=Jex_func(x);
global J_satex
global K_ex
y=J_satex*x/(x+K_ex); %See Eq. (5)
%y=x;

%To compare with Figure 2c of http://journal.frontiersin.org/article/10.3389/fnmol.2015.00006/full#B29
function r=alpha_minus_beta_cG(cG_left,Ca_left,cG_right,Ca_right,dx_left,dx_right);
%Do the interval [0,x0]
%Calculate currents for each spatiotemporal point
JcG=arrayfun(@a_b_nofeedback_func,cG_left);
%Jex=arrayfun(@Jex_func,Ca_left);
%Integrate over the length of the outer segment for each t
sum_left=trapz(transpose(JcG));
%sum_left=trapz(transpose(cG_left));
%Do the interval [x0,L]
%Calculate currents for each spatiotemporal point
JcG=arrayfun(@a_b_nofeedback_func,cG_right);
%Jex=arrayfun(@Jex_func,Ca_right);
%Integrate over the length of the outer segment for each t
sum_right=trapz(transpose(JcG));
%sum_right=trapz(transpose(cG_right));
tmp_r=(sum_left*dx_left+sum_right*dx_right);
%r=I_dark-2*(sum_left*dx_left);
%Take finite differences (time-derivative)
for i=1:length(tmp_r)-1
    dr(i)=(tmp_r(i+1)-tmp_r(i))*(200/(2.0*22));
end
r=dr;

function y=a_b_nofeedback_func(x);
global beta_dark
global alpha_dark
y=x;
%y=alpha_dark-beta_dark*x;
%global L_OS;
%y=-beta_dark*x;

%Modifying the example described in http://www.mathworks.com/help/matlab/ref/pdepe.html
%to solve the system of parabolic PDEs
function [u1,u2]=pde_call_right(x,t,option);
m = 0;
if option==0
    sol = pdepe(m,@pdex4pde_constantalpha,@pdex4ic,@pdex4bc_right,x,t);
else
    sol = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc_right,x,t);
end
u1 = sol(:,:,1);
u2 = sol(:,:,2);

function [u1,u2]=pde_call_left(x,t,option);
m = 0;
if option==0
    sol = pdepe(m,@pdex4pde_constantalpha,@pdex4ic,@pdex4bc_left,x,t);
else
    sol = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc_left,x,t);
end
u1 = sol(:,:,1);
u2 = sol(:,:,2);

function [u1,u2]=pde_call_alldisks(x,t,option);
m = 0;
if option==0
    sol = pdepe(m,@pdex4pde_alldisks_constantalpha,@pdex4ic,@pdex4bc_alldisks,x,t);
else
    sol = pdepe(m,@pdex4pde_alldisks,@pdex4ic,@pdex4bc_alldisks,x,t);
end
u1 = sol(:,:,1);
u2 = sol(:,:,2);

% --------------------------------------------------------------
function [c,f,s] = pdex4pde(x,t,u,DuDx);
global D_cG
global D_Ca
global beta_dark
global n_cyc
global K_cyc
global alpha_max
global cG_dark
global J_cGdark_f2
global FA2B
global J_satex
global K_ex
global n_cG
%global alpha_dark
%See Eqs. 1 and 2
c = [1; 1];
f = [D_cG; D_Ca] .* DuDx;
s = [alpha_max/(1+(u(2)/K_cyc)^n_cyc)-beta_dark*u(1); FA2B*(J_cGdark_f2*(u(1)/cG_dark)^n_cG-J_satex*u(2)/(u(2)+K_ex))];
%global Ca_dark
%s = [alpha_max/(1+(Ca_dark/K_cyc)^n_cyc)-beta_dark*u(1); FA2B*(J_cGdark_f2*(u(1)/cG_dark)^n_cG-J_satex*u(2)/(u(2)+K_ex))];

%Use Eq. (9) and incorporate a constant rate of cGMP synthesis
function [c,f,s] = pdex4pde_constantalpha(x,t,u,DuDx);
global D_cG
global D_Ca
global beta_dark
%global n_cyc
%global K_cyc
%global alpha_max
global cG_dark
global J_cGdark_f2
global FA2B
global J_satex
global K_ex
global alpha_dark
global n_cG
%See Eqs. 1 and 2
c = [1; 1];
f = [D_cG; D_Ca] .* DuDx;
s = [alpha_dark-beta_dark*u(1); FA2B*(J_cGdark_f2*(u(1)/cG_dark)^n_cG-J_satex*u(2)/(u(2)+K_ex))];

%No point sink for photon-induced cGMP hydrolysis.
%A term proportional to active PDE (E*(t)) is added to cGMP hydrolysis
%rate. If the initial spatial profile of the concentrations is constant,
%then the solutions should not vary spatially.
function [c,f,s] = pdex4pde_alldisks_constantalpha(x,t,u,DuDx);
global D_cG
global D_Ca
global beta_dark
%global n_cyc
%global K_cyc
%global alpha_max
global cG_dark
global J_cGdark_f2
global FA2B
global J_satex
global K_ex
global alpha_dark

global beta_idv
global R0
global k_R
global k_R2
global k_Rfraction
global k_E
global k_E2
global k_Efraction
global nu_RE

global n_cG

%See Eqs. 1 and 2
c = [1; 1];
f = [D_cG; D_Ca] .* DuDx;
expkRt=exp(-k_R*t);
expkR2t=exp(-k_R2*t);
expkEt=exp(-k_E*t);
expkE2t=exp(-k_E2*t);
E1=k_Rfraction*R0*nu_RE/(k_E-k_R)*(expkRt-expkEt)+(1-k_Rfraction)*R0*nu_RE/(k_E-k_R2)*(expkR2t-expkEt);
E2=k_Rfraction*R0*nu_RE/(k_E2-k_R)*(expkRt-expkE2t)+(1-k_Rfraction)*R0*nu_RE/(k_E2-k_R2)*(expkR2t-expkE2t);
s = [alpha_dark-(beta_dark + 0.5*beta_idv*(k_Efraction*E1+(1-k_Efraction)*E2) )*u(1); FA2B*(J_cGdark_f2*(u(1)/cG_dark)^n_cG-J_satex*u(2)/(u(2)+K_ex))];
%s = [alpha_dark-(beta_dark + 0.5*beta_idv*(k_Rfraction*R0*nu_RE/(k_E-k_R)*(exp(-k_R*t)-exp(-k_E*t))+(1-k_Rfraction)*R0*nu_RE/(k_E-k_R2)*(exp(-k_R2*t)-exp(-k_E*t))) )*u(1); FA2B*(J_cGdark_f2*(u(1)/cG_dark)^3-J_satex*u(2)/(u(2)+K_ex))];


function [c,f,s] = pdex4pde_alldisks(x,t,u,DuDx);
global D_cG
global D_Ca
global beta_dark
global n_cyc
global K_cyc
global alpha_max
global cG_dark
global J_cGdark_f2
global FA2B
global J_satex
global K_ex
%global alpha_dark

global beta_idv
global R0
global k_R
global k_R2
global k_Rfraction
global k_E
global k_E2
global k_Efraction
global nu_RE

global n_cG
global Ca_dark

%See Eqs. 1 and 2
c = [1; 1];
f = [D_cG; D_Ca] .* DuDx;
%f = [0; 0] .* DuDx; %Because diffusion does not matter (no spatial
%variation in x, setting the diffusion coefficients to zero gives the same
%result
expkRt=exp(-k_R*t);
expkR2t=exp(-k_R2*t);
expkEt=exp(-k_E*t);
expkE2t=exp(-k_E2*t);
if 1
    E1=k_Rfraction*R0*nu_RE/(k_E-k_R)*(expkRt-expkEt)+(1-k_Rfraction)*R0*nu_RE/(k_E-k_R2)*(expkR2t-expkEt);
else
    %Special combination
    E1=k_R*(R0*nu_RE/(k_E-k_R)*(expkRt-expkEt)-R0*nu_RE/(k_E-k_R2)*(expkR2t-expkEt))/(k_R2-k_R) + R0*nu_RE/(k_E-k_R)*(expkRt-expkEt);
    %To make sure
    %k_Efraction=1;
end
E2=k_Rfraction*R0*nu_RE/(k_E2-k_R)*(expkRt-expkE2t)+(1-k_Rfraction)*R0*nu_RE/(k_E2-k_R2)*(expkR2t-expkE2t);
%11/12/13 Put a (floor on [Ca2+]) ceiling on alpha in order to prolong Tsat for high intensities
%ff=1.0;
%if u(2)<ff*K_cyc
%    s = [alpha_max/(1+(ff)^n_cyc)-(beta_dark + 0.5*beta_idv*(k_Efraction*E1+(1-k_Efraction)*E2) )*u(1); FA2B*(J_cGdark_f2*(u(1)/cG_dark)^n_cG-J_satex*u(2)/(u(2)+K_ex))];
%else
    s = [alpha_max/(1+(u(2)/K_cyc)^n_cyc)-(beta_dark + 0.5*beta_idv*(k_Efraction*E1+(1-k_Efraction)*E2) )*u(1); FA2B*(J_cGdark_f2*(u(1)/cG_dark)^n_cG-J_satex*u(2)/(u(2)+K_ex))];
%end
%s = [alpha_max/(1+(u(2)/K_cyc)^n_cyc)-(beta_dark + 0.5*beta_idv*(k_Efraction*E1+(1-k_Efraction)*E2) )*u(1); FA2B*(J_cGdark_f2*(u(1)/cG_dark)^n_cG-J_satex*u(2)/(u(2)+K_ex))];
%s = [alpha_max/(1+(u(2)/K_cyc)^n_cyc)-(beta_dark + 0.5*beta_idv*(k_Rfraction*R0*nu_RE/(k_E-k_R)*(exp(-k_R*t)-exp(-k_E*t))+(1-k_Rfraction)*R0*nu_RE/(k_E-k_R2)*(exp(-k_R2*t)-exp(-k_E*t))) )*u(1); FA2B*(J_cGdark_f2*(u(1)/cG_dark)^3-J_satex*u(2)/(u(2)+K_ex))];

% --------------------------------------------------------------
function u0 = pdex4ic(x);
global cG_dark
global Ca_dark
u0 = [cG_dark; Ca_dark];
% --------------------------------------------------------------
%Use this BC for the domain to the right (x>x0) of the isomerization point.
function [pl,ql,pr,qr] = pdex4bc_right(xl,ul,xr,ur,t);
global delta_beta4
global R0
global k_R
global k_R2
global k_Rfraction
global k_E
global k_E2
global k_Efraction
global nu_RE
expkRt=exp(-k_R*t);
expkR2t=exp(-k_R2*t);
expkEt=exp(-k_E*t);
expkE2t=exp(-k_E2*t);
E1=k_Rfraction*R0*nu_RE/(k_E-k_R)*(expkRt-expkEt)+(1-k_Rfraction)*R0*nu_RE/(k_E-k_R2)*(expkR2t-expkEt);
E2=k_Rfraction*R0*nu_RE/(k_E2-k_R)*(expkRt-expkE2t)+(1-k_Rfraction)*R0*nu_RE/(k_E2-k_R2)*(expkR2t-expkE2t);
pl=[-delta_beta4*(k_Efraction*E1+(1-k_Efraction)*E2)*ul(1); 0];
%disp testul1
%ul(1) %TODO Calculate derivative from solution and compare with this BC
%pl = [-delta_beta4*(k_Rfraction*R0*nu_RE/(k_E-k_R)*(exp(-k_R*t)-exp(-k_E*t))+(1-k_Rfraction)*R0*nu_RE/(k_E-k_R2)*(exp(-k_R2*t)-exp(-k_E*t)))*ul(1); 0];
%pl = [-delta_beta4*R0*nu_RE/(k_E-k_R)*(exp(-k_R*t)-exp(-k_E*t))*ul(1); 0]; %See Eq. (6). No diffusion constant D_cG here, see http://www.mathworks.com/help/matlab/ref/pdepe.html
%pl = [0; 0];
ql = [1; 1];
pr = [0; 0]; 
qr = [1; 1]; 

%Use this BC for the domain to the left (x<x0) of the isomerization point.
function [pl,ql,pr,qr] = pdex4bc_left(xl,ul,xr,ur,t);
global delta_beta4
global R0
global k_R
global k_R2
global k_Rfraction
global k_E
global k_E2
global k_Efraction
global nu_RE
expkRt=exp(-k_R*t);
expkR2t=exp(-k_R2*t);
expkEt=exp(-k_E*t);
expkE2t=exp(-k_E2*t);
E1=k_Rfraction*R0*nu_RE/(k_E-k_R)*(expkRt-expkEt)+(1-k_Rfraction)*R0*nu_RE/(k_E-k_R2)*(expkR2t-expkEt);
E2=k_Rfraction*R0*nu_RE/(k_E2-k_R)*(expkRt-expkE2t)+(1-k_Rfraction)*R0*nu_RE/(k_E2-k_R2)*(expkR2t-expkE2t);
pr=[delta_beta4*(k_Efraction*E1+(1-k_Efraction)*E2)*ur(1); 0];
%disp testur1
%ur(1) %TODO Calculate derivative from solution and compare with this BC
%pr = [delta_beta4*(k_Rfraction*R0*nu_RE/(k_E-k_R)*(exp(-k_R*t)-exp(-k_E*t))+(1-k_Rfraction)*R0*nu_RE/(k_E-k_R2)*(exp(-k_R2*t)-exp(-k_E*t)))*ur(1); 0];
%pr = [delta_beta4*R0*nu_RE/(k_E-k_R)*(exp(-k_R*t)-exp(-k_E*t))*ur(1); 0]; %See Eq. (6)
%pl = [0; 0];
qr = [1; 1];
pl = [0; 0]; 
ql = [1; 1]; 

%Set the fluxes to zero at the boundaries
function [pl,ql,pr,qr] = pdex4bc_alldisks(xl,ul,xr,ur,t);
pr = [0; 0];
qr = [1; 1]; 
pl = [0; 0]; 
ql = [1; 1]; 

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_k_E_Callback(hObject, eventdata, handles)
% hObject    handle to edit_k_E (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_k_E as text
%        str2double(get(hObject,'String')) returns contents of edit_k_E as a double


% --- Executes during object creation, after setting all properties.
function edit_k_E_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_k_E (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_k_R_Callback(hObject, eventdata, handles)
% hObject    handle to edit_k_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_k_R as text
%        str2double(get(hObject,'String')) returns contents of edit_k_R as a double


% --- Executes during object creation, after setting all properties.
function edit_k_R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_k_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_nu_RE_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nu_RE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nu_RE as text
%        str2double(get(hObject,'String')) returns contents of edit_nu_RE as a double


% --- Executes during object creation, after setting all properties.
function edit_nu_RE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nu_RE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_timediv_Callback(hObject, eventdata, handles)
% hObject    handle to edit_timediv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_timediv as text
%        str2double(get(hObject,'String')) returns contents of edit_timediv as a double


% --- Executes during object creation, after setting all properties.
function edit_timediv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_timediv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_x0leftdiv_Callback(hObject, eventdata, handles)
% hObject    handle to edit_x0leftdiv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_x0leftdiv as text
%        str2double(get(hObject,'String')) returns contents of edit_x0leftdiv as a double


% --- Executes during object creation, after setting all properties.
function edit_x0leftdiv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_x0leftdiv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_x0rightdiv_Callback(hObject, eventdata, handles)
% hObject    handle to edit_x0rightdiv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_x0rightdiv as text
%        str2double(get(hObject,'String')) returns contents of edit_x0rightdiv as a double


% --- Executes during object creation, after setting all properties.
function edit_x0rightdiv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_x0rightdiv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_timepoint_cG_Callback(hObject, eventdata, handles)
% hObject    handle to edit_timepoint_cG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_timepoint_cG as text
%        str2double(get(hObject,'String')) returns contents of edit_timepoint_cG as a double


% --- Executes during object creation, after setting all properties.
function edit_timepoint_cG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_timepoint_cG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_timeafter_Callback(hObject, eventdata, handles)
% hObject    handle to edit_timeafter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_timeafter as text
%        str2double(get(hObject,'String')) returns contents of edit_timeafter as a double


% --- Executes during object creation, after setting all properties.
function edit_timeafter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_timeafter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_alpha_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_alpha_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_alpha_max as text
%        str2double(get(hObject,'String')) returns contents of edit_alpha_max as a double


% --- Executes during object creation, after setting all properties.
function edit_alpha_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_alpha_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_alpha_dark_Callback(hObject, eventdata, handles)
% hObject    handle to edit_alpha_dark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_alpha_dark as text
%        str2double(get(hObject,'String')) returns contents of edit_alpha_dark as a double


% --- Executes during object creation, after setting all properties.
function edit_alpha_dark_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_alpha_dark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_n_cyc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_n_cyc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_n_cyc as text
%        str2double(get(hObject,'String')) returns contents of edit_n_cyc as a double


% --- Executes during object creation, after setting all properties.
function edit_n_cyc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_n_cyc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_K_cyc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_K_cyc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_K_cyc as text
%        str2double(get(hObject,'String')) returns contents of edit_K_cyc as a double


% --- Executes during object creation, after setting all properties.
function edit_K_cyc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_K_cyc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Rstar_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Rstar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Rstar as text
%        str2double(get(hObject,'String')) returns contents of edit_Rstar as a double


% --- Executes during object creation, after setting all properties.
function edit_Rstar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Rstar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_B_Ca_Callback(hObject, eventdata, handles)
% hObject    handle to edit_B_Ca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_B_Ca as text
%        str2double(get(hObject,'String')) returns contents of edit_B_Ca as a double


% --- Executes during object creation, after setting all properties.
function edit_B_Ca_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_B_Ca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_f_Ca_Callback(hObject, eventdata, handles)
% hObject    handle to edit_f_Ca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_f_Ca as text
%        str2double(get(hObject,'String')) returns contents of edit_f_Ca as a double


% --- Executes during object creation, after setting all properties.
function edit_f_Ca_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_f_Ca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_cG_dark_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cG_dark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cG_dark as text
%        str2double(get(hObject,'String')) returns contents of edit_cG_dark as a double


% --- Executes during object creation, after setting all properties.
function edit_cG_dark_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cG_dark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Ca_dark_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Ca_dark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Ca_dark as text
%        str2double(get(hObject,'String')) returns contents of edit_Ca_dark as a double


% --- Executes during object creation, after setting all properties.
function edit_Ca_dark_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Ca_dark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_beta_idv_Callback(hObject, eventdata, handles)
% hObject    handle to edit_beta_idv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_beta_idv as text
%        str2double(get(hObject,'String')) returns contents of edit_beta_idv as a double


% --- Executes during object creation, after setting all properties.
function edit_beta_idv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_beta_idv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_beta_dark_Callback(hObject, eventdata, handles)
% hObject    handle to edit_beta_dark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_beta_dark as text
%        str2double(get(hObject,'String')) returns contents of edit_beta_dark as a double


% --- Executes during object creation, after setting all properties.
function edit_beta_dark_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_beta_dark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_J_cGdark_Callback(hObject, eventdata, handles)
% hObject    handle to edit_J_cGdark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_J_cGdark as text
%        str2double(get(hObject,'String')) returns contents of edit_J_cGdark as a double


% --- Executes during object creation, after setting all properties.
function edit_J_cGdark_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_J_cGdark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_J_satex_Callback(hObject, eventdata, handles)
% hObject    handle to edit_J_satex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_J_satex as text
%        str2double(get(hObject,'String')) returns contents of edit_J_satex as a double


% --- Executes during object creation, after setting all properties.
function edit_J_satex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_J_satex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_K_ex_Callback(hObject, eventdata, handles)
% hObject    handle to edit_K_ex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_K_ex as text
%        str2double(get(hObject,'String')) returns contents of edit_K_ex as a double


% --- Executes during object creation, after setting all properties.
function edit_K_ex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_K_ex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_I_dark_Callback(hObject, eventdata, handles)
% hObject    handle to edit_I_dark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_I_dark as text
%        str2double(get(hObject,'String')) returns contents of edit_I_dark as a double


% --- Executes during object creation, after setting all properties.
function edit_I_dark_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_I_dark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_alldisks.
function checkbox_alldisks_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_alldisks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_alldisks



function edit_kR2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_kR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_kR2 as text
%        str2double(get(hObject,'String')) returns contents of edit_kR2 as a double


% --- Executes during object creation, after setting all properties.
function edit_kR2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_kR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_kR1fraction_Callback(hObject, eventdata, handles)
% hObject    handle to edit_kR1fraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_kR1fraction as text
%        str2double(get(hObject,'String')) returns contents of edit_kR1fraction as a double


% --- Executes during object creation, after setting all properties.
function edit_kR1fraction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_kR1fraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_k_E2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_k_E2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_k_E2 as text
%        str2double(get(hObject,'String')) returns contents of edit_k_E2 as a double


% --- Executes during object creation, after setting all properties.
function edit_k_E2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_k_E2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_kE1fraction_Callback(hObject, eventdata, handles)
% hObject    handle to edit_kE1fraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_kE1fraction as text
%        str2double(get(hObject,'String')) returns contents of edit_kE1fraction as a double


% --- Executes during object creation, after setting all properties.
function edit_kE1fraction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_kE1fraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
