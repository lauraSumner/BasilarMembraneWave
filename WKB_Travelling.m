% Calculate the amplitude of BM waves based on the WKB approximation and information obtained by the experimental work of Emadi et al.

% Code shows the effects at different stimili frequency.

%close all
clear all

%The following variables set the parameters which define the BM wave. These
%values apply to the cochlea of a small rodent (guinea pig/chinchilla) and
%were fixed according to experimental data from such species. For a human
%cochlea, the values should be altered accordingly.
% Parameters: height of channel, density of fluid, mass at x=0, stiffness
% at x=0, stiffness exp decay factor, mass exp increase factor.
h=0.5e-3;   %Height of channel
rho=1000;   %Fluid density
K0=10;      %Stiffness function amplitude
m0=32e-11;  %Mass function amplitude
a=501;      %Exp decay factor
b=101;      %Exp growth factor
L=0.011;    %Length of cochlea
eta=1e-6;   %Fluid viscosity

dampFac=[10e-8, 5e-7, 2e-6 ]; %This has been fixed accordingly to experiemtal data to best match the amplitudes
%of wave envelopes between SPL (60, 80 100). The tuning is however stll quite sharp.
arb=8e5;
pressure=[0.02,  0.2, 2];%pressure scale factors between SPLs.
%Geometric parameters D and W used to calcualte the area, A  of an 8um section
%of the membrane.
D=168e-6;
W=204e-6;
A=(W+D)*(8e-6)/2;
BMpoints=979; %chose to match the number of BM points, found in the boundary file (number of faces/2  + 1)
% Time parameters should be chosen to match the simulation
T=0.002;
dt=5e-6;
t=[0:dt:T];
timesteps=length(t);
x=linspace(0, L,BMpoints);
N=length(x);
H=L/N;
freq=[10000]; %can change this to an array across different freqencies for quicker output.
omega=2*pi*freq;


% Empty arrays
arrayFunc=zeros(1, N);
final=zeros(length(omega), N);
newmass= zeros(length(omega), N);
kVector= zeros(length(omega), N);

% exponential functions
K=K0*exp(-a*x);
m=m0*exp(b*x);

%loop over SPLs
for count= 1:3
    
    arbitrary=arb*pressure(count);
    damp=dampFac(count);
    
    for i= 1:length(omega)
        omega1=omega(i);
        om=num2str(omega1);
        
        
        C=sqrt(-2*1i*rho*omega1/h);
        cArray(i)=C;
        
        %Impedence
        Z1= (-1i*K/omega1).*(1-((omega1^2)*m./K))/A + damp/A;
        
        
        %Prefactor
        V=arbitrary*Z1.^(-3/4)./C;
        vArray(i,:)= V;
        
        %Wave Vector
        kVector= C./sqrt(Z1);
        
        %Integration function
        
        for k= 2:N
            
            func= trapz(x(1:k), kVector(1:k));
            
            arrayFunc(k)= func;
            
        end
        
        kVectorArray(i, :) = kVector;
        
        %Full Spatial Function
        final(i, :)= V.*exp(-1i*arrayFunc);
        displacement(i,:)=final(i,:)./omega1;
        velEnv(i,:,count)=abs(final(i,:));
        dispEnv(i,:,count)=abs(displacement(i,:));
        
    end
    %Full Time-depedent function.
    for i= 1:timesteps
        time=t(i);
        y= displacement(1,:) .* (exp(1i*omega1*time));
        
        yMatrix(i,:) = real(y);
        
        x1=zeros(N); % need the x file to be a zero file
        
    end
    yMatrix=flip(yMatrix,2);
    
    dimMatrix(1,1)=length(t);
    dimMatrix(1,2)=length(x);
    
    
    
    save('dimensions_y.txt','dimMatrix','-ascii');
    save('matrix_y.txt','yMatrix','-ascii');
    save('matrix_x.txt','x1','-ascii');
    
    
    
end

