%Author:Amit Herman, third year student 2017

function [ ] = Q1()
fprintf('Welcome\n');
fprintf('Please select one of the options:\n');
fprintf('1 for Bet\n');
fprintf('2 for Gimel\n');
fprintf('3 for Dalet\n');
fprintf('4 for He\n');
n = input('Enter a number: ');
switch n
    case 1
    ex = 1;
    ey = 0;
    ez = 0;
    L1 = 11/8;
    R = 1;
    Nz = 50;
    Nphi = 60;
    Bet(Nz, Nphi, R, ex, ey ,ez, L1);
    
    case 2
    ex = 0;
    ey = 0;
    ez = 1;
    L1 = 1;
    R = 1;
    Nz = 50;
    Nphi = 60;
    Gimel(Nz, Nphi, R, ex, ey ,ez, L1);
    case 3
    ex = 1;
    ey = 0;
    ez = 1;
    L1 = 1;
    R = 1;
    Nz = 50;
    Nphi = 60;
    Dalet(Nz, Nphi,R,ex, ey, ez, L1, 10);
    
        case 4
    ex = 1;
    ey = 0;
    ez = 1;
    L1 = 1;
    R = 1;
    Nz = 50;
    Nphi = 60;
    He(Nz, Nphi,R,ex, ey, ez, L1, 10);
end
end





function [ Sigma ] = FindSigma(Nz, Nphi, r, ex, ey, ez, L1)
%Input: Division order in phi direction and z direction,radius , external field  in each direction and L1
%the height of the shell
%Output: The function returns charge distribution in each element in the
%shell


N_1 = Nz*Nphi;
epsilon_0 = 8.8541 * 10^-12;
delta_phi = 2*pi/Nphi;
delta_z = L1/Nz;
b_phi = delta_phi*r;
b_z = delta_z;
alpha = b_phi/b_z;


%create series of vectors
m = 1:Nz;
n = 1:Nphi;
[M,N] = meshgrid(m, n);
Vecs_x = r*cos(delta_phi.*(N - 0.5));
Vecs_y =  r*sin(delta_phi.*(N - 0.5));
Vecs_z = -1*L1/2 + delta_z.*(M-0.5);

%create matrice L now
A = zeros(N_1);
for i = 1:N_1 % m
    for j = 1:N_1 % n
       % iter = iter + 1;
        if i~=j           
        A(i, j) =  (1/(4*pi*epsilon_0)) *b_z*b_phi/norm([Vecs_x(j) Vecs_y(j) Vecs_z(j)] -[Vecs_x(i) Vecs_y(i) Vecs_z(i)]) ;
        else 
            A(i, j) = (1/(4*pi*epsilon_0)) *log((sqrt(alpha^2 +1) + alpha)/(sqrt(alpha^2 +1) - alpha))*b_z + ...
                      (1/(4*pi*epsilon_0)) *log((sqrt(alpha^2 +1) + 1)/(sqrt(alpha^2 +1) - 1))*b_phi;
        end
    end
end

%The following lines handle the n+1'th equation!
A(N_1+1, 1:N_1) = b_phi*b_z;
A(N_1+1, N_1+1) = 0;
A(1:N_1, N_1+1) = -1;
%Now we will make vector v
V = zeros(N_1, 1);

for i = 1:N_1
    V(i) = [Vecs_x(i) Vecs_y(i) Vecs_z(i)]*transpose([ex ey ez]);
end 
    V(N_1 + 1) = 0;

%And now at least :)
Sigma = A\V;
       
end

function [Dist ] = Find_distribution_z(Sigma, z0, Nz, Nphi, L1)
%Input:Charge distribution on the shell, z value , division order on z,phi and L1 the height.
%output:The function returns the distribution on the closest curve when z=z0.
b_z = L1/Nz;

N = floor( (z0+ L1/2)/b_z + 0.5);% the nearest index 
if(N==0)
    N = 1;
end
Dist = Sigma(1 + (N-1)*Nphi :N*Nphi);

end

function [Dist ] = Find_distribution_phi(Sigma, phi0, Nz, Nphi,r )
%Input:Charge distribution on the shell, phi value , division order on z,phi and radius.
%output:The function returns the distribution on the closest line to phi = phi0.
b_phi = (2*pi/Nphi)*r;

N = floor( ((phi0/b_phi)+ 0.5));% the nearest index 
if(N==0)
    N = 1;
end
Dist = Sigma(N:Nphi:Nphi*Nz);

end


function [Px] = Find_Dipole_momx(Sigma, Nz, Nphi, L1, radius)
%Input:Charge distribution on the shell, division order of z,of phi, L1 the height and radius.
%Output:Dipole moment in X direction

delta_phi = 2*pi/Nphi;
b_z = L1/Nz;
b_phi = (2*pi/Nphi)*radius;
ds = b_z*b_phi;% Area element
Px = 0;


for i = 1:Nz*Nphi
    n = mod(i-1, Nphi)+ 1;
    phi = delta_phi*(n - 0.5);
    Px  = Px + radius*cos(phi)*Sigma(i)*ds;
end

end

function [Pz] = Find_Dipole_momz(Sigma, Nz, Nphi, L1, radius)
%Input:Charge distribution on the shell, division order of z,of phi, L1 the height and radius.
%Output:Dipole moment in Z direction

b_z = L1/Nz; % also delta z
b_phi = (2*pi/Nphi)*radius;
ds = b_z*b_phi;% Area element
Pz = 0;

for i = 1:Nz*Nphi
    n = floor(i-1/Nphi) + 1;
    z = -L1/2 + b_z*(n-0.5);
    Pz  = Pz + z*Sigma(i)*ds;
end

end

function [alpha] = calc_polar(dipole, ex_f)
%Input:Dipole moment, external electric field in the same direction.
%output:The function calculates the polarizability in that direction

epsilon_0 = 8.8541 * 10^-12;
alpha = dipole/(epsilon_0*ex_f);
end


function [] = Bet(Nz, Nphi, r, ex, ey, ez, L1)
%This function solves Bet
Sig = FindSigma(Nz, Nphi, r, ex, ey, ez, L1);
fprintf('Excerice Bet\n');
fprintf('Please select one of the options:\n');
fprintf('1 for Bet(2)\n');
fprintf('2 for Bet(3)\n');
fprintf('999 to stop\n');
n = input('Enter a number: ');
while(n~=999)
    
switch n
    case 1
        Bet_2(Sig, Nphi, Nz, L1, r);
    case 2
        Bet_3(Sig, Nz, Nphi, L1, r, ex);

end


fprintf('\nPlease select one of the options:\n');
fprintf('1 for Bet(2)\n');
fprintf('2 for Bet(3)\n');
fprintf('999 to stop\n');
n = input('Enter a number: ');
end



end

function [] = Bet_2(Sigma, Nphi, Nz, L1, R)
%This function recieves distribution, Nphi and Nz and solves Bet_2


phi = 2*pi/Nphi*((1:Nphi) - 0.5);
z = -L1/2 + L1/Nz*((1:Nz) - 0.5);
fprintf('You have selected Beth 2\n');
fprintf('Please select one of the options:\n');
fprintf('1 for phi = 0\n');
fprintf('2 for phi = pi/2\n');
fprintf('3 for phi = pi\n');
fprintf('4 for z = -L1/4\n');
fprintf('5 for z = 0\n');
fprintf('6 for z = L1/4\n');
fprintf('999 to stop\n');
n = input('Enter a number: ');

while(n~=999)
switch n
    case 1
       D  = Find_distribution_phi(Sigma, 0, Nz, Nphi, R);
       figure;
       plot(z, D);
       
      title('Phi = 0');
    xlabel('z(m*10^-5)');
    ylabel('צפיפות מטען משטחית');
    case 2
        D  = Find_distribution_phi(Sigma, pi/2, Nz, Nphi, R);
         figure;
          plot(z, D);
          title('Phi = Pi/2');
         xlabel('z(m*10^-5)');
       ylabel('צפיפות מטען משטחית');
    case 3
        D  = Find_distribution_phi(Sigma, pi, Nz, Nphi, R);
         figure;
       plot(z, D);
            title('Phi = Pi');
         xlabel('z(m*10^-5)');
       ylabel('צפיפות מטען משטחית');
    case 4
        D = Find_distribution_z(Sigma, -L1/4, Nz, Nphi, L1);
         figure;
        plot(phi, D);
                    title('Z = -L1/4');
         xlabel('z(m*10^-5)');
       ylabel('צפיפות מטען משטחית');
     case 5
        D = Find_distribution_z(Sigma, 0, Nz, Nphi, L1);
         figure;
        plot(phi, D);
                    title('Z = 0');
         xlabel('z(m*10^-5)');
       ylabel('צפיפות מטען משטחית');
    case 6
        D = Find_distribution_z(Sigma, L1/4, Nz, Nphi, L1);
         figure;
        plot(phi, D);
                   title('Z = L1/4');
         xlabel('z(m*10^-5)');
       ylabel('צפיפות מטען משטחית');
end
n = input('Enter a number: ');
end
fprintf('\n');

end


function[] = Bet_3(Sigma, Nz, Nphi, L1,R, ex)
%This function solves Bet 3.
Dipole = Find_Dipole_momx(Sigma, Nz, Nphi,L1, R);
alpha = calc_polar(Dipole, ex);%Polarizability
fprintf('The value of the dipole moment is %d [C*m*10^-5]\n', Dipole );
fprintf('The Polarizability is %d\n', alpha );

fprintf('\n');
end


function [] = Gimel(Nz, Nphi, r, ex, ey, ez, L1)
%This function solves Gimel
Sig = FindSigma(Nz, Nphi, r, ex, ey, ez, L1);
fprintf('Excerice Gimel\n');
fprintf('Please select one of the options:\n');
fprintf('1 for Gimel(2)\n');
fprintf('2 for Gimel(3)\n');
fprintf('999 to stop\n');
n = input('Enter a number: ');

while(n~=999)
switch n
    case 1
        Gimel_2(Sig, Nphi, Nz, L1, r);
    case 2
        Gimel_3(Sig, Nz, Nphi, L1, r, ez);

end


fprintf('Please select one of the options:\n');
fprintf('1 for Gimel(2)\n');
fprintf('2 for Gimel(3)\n');
fprintf('999 to stop\n');
n = input('Enter a number: ');
end

end

function [] = Gimel_2(Sigma, Nphi, Nz, L1, R)
%This function recieves distribution, Nphi and Nz and solves Gimel_2


phi = 2*pi/Nphi*((1:Nphi) - 0.5);
z = -L1/2 + L1/Nz*((1:Nz) - 0.5);
fprintf('\nYou have selected Gimel 2\n');
fprintf('Please select one of the options:\n');
fprintf('1 for phi = 0\n');
fprintf('2 for phi = pi/2\n');
fprintf('3 for phi = pi\n');
fprintf('4 for z = -L1/4\n');
fprintf('5 for z = 0\n');
fprintf('6 for z = L1/4\n');
fprintf('999 to stop\n');
n = input('Enter a number: ');

while(n~=999)
switch n
       case 1
       D  = Find_distribution_phi(Sigma, 0, Nz, Nphi, R);
       figure;
       plot(z, D);
       
      title('Phi = 0');
    xlabel('z(m*10^-5)');
    ylabel('צפיפות מטען משטחית');
    case 2
        D  = Find_distribution_phi(Sigma, pi/2, Nz, Nphi, R);
         figure;
          plot(z, D);
          title('Phi = Pi/2');
         xlabel('z(m*10^-5)');
       ylabel('צפיפות מטען משטחית');
    case 3
        D  = Find_distribution_phi(Sigma, pi, Nz, Nphi, R);
         figure;
       plot(z, D);
            title('Phi = Pi');
         xlabel('z(m*10^-5)');
       ylabel('צפיפות מטען משטחית');
    case 4
        D = Find_distribution_z(Sigma, -L1/4, Nz, Nphi, L1);
         figure;
        plot(phi, D);
        ylim([-50*10^-11 50*10^-11]);
                    title('Z = -L1/4');
         xlabel('z(m*10^-5)');
       ylabel('צפיפות מטען משטחית');
     case 5
        D = Find_distribution_z(Sigma, 0, Nz, Nphi, L1);
         figure;
        plot(phi, D);
        ylim([-50*10^-11 50*10^-11]);
                    title('Z = 0');
         xlabel('z(m*10^-5)');
       ylabel('צפיפות מטען משטחית');
    case 6
        D = Find_distribution_z(Sigma, L1/4, Nz, Nphi, L1);
         figure;
        plot(phi, D);
       ylim([-50*10^-11 50*10^-11]);

                   title('Z = L1/4');
         xlabel('z(m*10^-5)');
       ylabel('צפיפות מטען משטחית');

end
n = input('Enter a number: ');
end
fprintf('\n');
end


function[] = Gimel_3(Sigma, Nz, Nphi, L1,R, ex)
%This function solves Gimel 3.
Dipole = Find_Dipole_momz(Sigma, Nz, Nphi,L1, R);
alpha = calc_polar(Dipole, ex);%Polarizability
fprintf('The value of the dipole moment is %d [C*m*10^-5]\n', Dipole );
fprintf('The Polarizability is %d \n', alpha );
fprintf('\n');
end


function[] = Dalet(Nz, Nphi, r, ex, ey, ez, L1A, n)
%The function recieves the division orders in z and phi, radius, external
%forces, L1A and n which is used to divide L1 axis in this excercise
%This function solves Dalet
%This may take some time so be patient please.
i=0;
alpha_xx = zeros(1, n+1);
alpha_zz = zeros(1,n+1);
delta_L1 = (2*L1A - 0.5*L1A)/n;
for i = 0:n
    Sig = FindSigma(Nz, Nphi, r, ex, ey, ez,L1A/2 + delta_L1*i);
    Px = Find_Dipole_momx(Sig, Nz, Nphi, L1A/2 + delta_L1*i, r);
    Pz = Find_Dipole_momz(Sig, Nz, Nphi, L1A/2 + delta_L1*i, r);
    
    alpha_xx(i+1) = calc_polar(Px, ex);
    alpha_zz(i+1) = calc_polar(Pz, ez);
end

L = L1A/2 + delta_L1*(0:n);
fprintf('What would you like to do?\n');
fprintf('1 for Alpha xx\n');
fprintf('2 for Alpha zz\n');
fprintf('999 to stop\n');
n = input('Enter a number: ');

while(n~=999)
switch n
    case 1
       figure;
       plot(L, alpha_xx);
    case 2
       figure;
       plot(L, alpha_zz);
end
fprintf('\nWhat would you like to do?\n');
fprintf('1 for Alpha xx\n');
fprintf('2 for Alpha zz\n');
fprintf('999 to stop\n');
n = input('Enter a number: ');
end


end



function[] = He(Nz, Nphi, r, ex, ey, ez, L1, n)
%The function recieves the division orders in z and phi, radius, external
%forces, L1 and n which is used to divide r axis in this excercise
%This function solves He
%This may take some time so be patient please.
i=0;
alpha_xx = zeros(1, n+1);
alpha_zz = zeros(1,n+1);
delta_R = (2*r - 0.5*r)/n;
for i = 0:n
    Sig = FindSigma(Nz, Nphi, 0.5*r + delta_R*i, ex, ey, ez,L1);
    Px = Find_Dipole_momx(Sig, Nz, Nphi, L1, 0.5*r + delta_R*i);
    Pz = Find_Dipole_momz(Sig, Nz, Nphi, L1, 0.5*r + delta_R*i);
    
    alpha_xx(i+1) = calc_polar(Px, ex);
    alpha_zz(i+1) = calc_polar(Pz, ez);
end

R = r/2 + delta_R*(0:n);
fprintf('What would you like to do?\n');
fprintf('1 for Alpha xx\n');
fprintf('2 for Alpha zz\n');
fprintf('999 to stop\n');
n = input('Enter a number: ');

while(n~=999)
switch n
    case 1
       figure;
       plot(R, alpha_xx);
                          title('Polarizabiltity - x');
         xlabel('Radius');
       ylabel('Polarizability');
    case 2
       figure;
       plot(R, alpha_zz);
                                 title('Polarizabiltity - z');
         xlabel('Radius');
       ylabel('Polarizability');
end
fprintf('\nWhat would you like to do?\n');
fprintf('1 for Alpha xx\n');
fprintf('2 for Alpha zz\n');
fprintf('999 to stop\n');
n = input('Enter a number: ');
end


end

