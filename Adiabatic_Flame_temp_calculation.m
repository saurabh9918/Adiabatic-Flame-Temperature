%% This is the Program to calculate the Adiabatic Flame temperature considering the Constant Pressure process
%% Units are to be taken according to the conditions:
%% Temperature in Kelvin  Unit of enthalpy is Joule and accordingly for other take according to indicated
%% Section 1- Taking the Input
x =input('Enter the value of integer value of x for deciding the carbon');
y =input('Enter the value of integer value of y for deciding the Hydrogen');
CpCxHy =input('Enter the value of Specfic Heat for Hydorcarbon considering it to be constant');
T1= input('Enter the sensible Heating temperature in Kelvin for the fuel and air');
dHfCxHy =input('Enter the Enthalpy of formation of Hydrocarbon at 298 K');
n= input('Enter the number of samples for graph plotting');
for i = 1:n
    phi =input('Enter the Equivalence ratio'); %taking equivalence ratio accordingly the mixture is lean or rich
    % consider the Hydrocarbon as CxHy we will take the value of x and y
    phiA(i)= phi;
%% Section 2- Assigning some values and forming the Algebaric Equation for Adiabatic Flame temperature
% Assuming Cp= A+BT+(C/T)
   dHfO2= 1;
   dHfN2= 2;
   dHfCO2= 3;
   dHfH2O= 4;
   C1= dHfCxHy+dHfO2+dHfN2;
   C2= dHfCO2+dHfH2O+dHfO2+dHfN2;
   C= C1-C2;
   n1= 1;
   n2= (x+(y/4))/phi;
   n3= (3.76/phi)*(x+(y/4));
   n4= x;
   n5= y/2;
   n6= n3;
   n7= ((1/phi)-1)*(x+(y/4));
   lamda1= C+n1*CpCxHy*(T1-298);
   lamda2= T1-298;
   lamda3= (T1*T1)-(298*298);
   lamda4= (1/T1)-(1/298);
%Now assigning all the coefficients of Cp for each element in a array, it
%can be done sepreatley
   CpO2= [4.184*0 4.184*(7.16e-3) -4.184*(0.4e5)];
   CpN2= [4.184*6.66 4.184*(1.02e-3) -4.184*(0)];
   CpCO2= [4.184*10.55 4.184*(2.16e-3) -4.184*(0.8e5)];
   CpH2O= [4.184*7.17 4.184*(2.56e-3) -4.184*(2.04e5)];
   p= lamda1+n2*((CpO2(1)*lamda2)+(CpO2(2)*lamda3)+(CpO2(3)*lamda4))+n3*((CpN2(1)*lamda2)+(CpN2(2)*lamda3)+(CpN2(3)*lamda4));
   px= n4*CpCO2(1)+n5*CpH2O(1)+n6*CpN2(1)+n7*CpO2(1);
   py= n4*CpCO2(2)+n5*CpH2O(2)+n6*CpN2(2)+n7*CpO2(2);
   pz= n4*CpCO2(3)+n5*CpH2O(3)+n6*CpN2(3)+n7*CpO2(3);
   sai= p+(298*px)+(298*298)*(py/2)+(pz/298);
%% Section3- Initializing the Array for the Cubic Equation
   cubice= [py/2 px -sai -pz];
%% Section4- Solving using Newton raphson method
   dcubice= [1.5*py 2*px -sai];
   guess= input('Enter the guess value in kelvin'); %getting the Initial guess value
   relax= input('Enter the relaxation factor'); %getting the relaxation factor
   fguess= cubice(1)*guess*guess*guess+cubice(2)*guess*guess+cubice(3)*guess+cubice(4);
   dfguess= dcubice(1)*guess*guess+dcubice(2)*guess+dcubice(3);
   x1= guess-(relax)*(fguess/dfguess);
   e =input('Enter the accuracy required in fraction');
   e1= abs(x1-guess);
   while e1>=e
    fguess= cubice(1)*x1*x1*x1+cubice(2)*x1*x1+cubice(3)*x1+cubice(4);
    dfguess= dcubice(1)*x1*x1+dcubice(2)*x1+dcubice(3);
    x2= x1-(relax)*(fguess/dfguess);
    e1= abs(x1-x2);
    x1= x2;
   end
   disp(x1);
   disp(e1);
   TaA(i)= x1;
%% The Program ends and the final error is noted whether the convergence is attained or not
end
plot(TaA,phiA); %graph for the releation



    


