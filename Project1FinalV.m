% Project I: Using Iteration Methods to Understand Fractal Geometry
% Abby Saenz, Cecilia Rodarte, Dannie Kiel, Taylor Reeder

clc; clear; close all;

global cval  
cval = [-1.25, 0.36 + 0.1*(1i), -.123 - .745*(1i), -0.70176 - 0.2842*(1i)];

%% Part I - An Introduction to Fractals
clearvars -except cval

phi = @(z) z^2;
x = linspace(-1,1,200);
y = linspace(-1,1,200);
M = ones(200,200);

for i = 1:200
    for j = 1:200
        clear z;
        z = x(i) + 1i*y(j);
        for k = 1:100
            z(k+1) = phi(z(k));
            if abs(z(k+1)) > 2
                M(i,j) = 2;
            end
        end
    end
end

% figure
% colormap([0 0 1; 1 1 1]);
% image( [-1 1], [-1 1], M)
% axis('equal')
% axis([ -1.5 1.5 -1.5 1.5])

%% Part II - Generating Other Examples for Various Complex Values c
clearvars -except cval

phi = @(z,cval) z^2 + cval;
x = linspace(-1.8,1.8,500);
y = linspace(-1.25,1.25,500);
M = cell(length(cval),1);

for k = 1:length(cval)
    M{k} = ones(length(x),length(y));
    for r = 1:length(x)
        for i = 1:length(y)
            clear z;
            z = x(r) + 1i*y(i);
            for j = 1:100
                z(j+1) = phi(z(j),cval(k));
                if abs(z(j+1)) > 2
                    M{k}(r,i) = 2;
                    break;
                end
            end
        end
    end
end

% for i = 1:length(cval)
%     M{i} = M{i}';
%     figure
%     hold on
%     colormap([0 0 1; 1 1 1]);
%     image( [-1 1], [-1 1], M{i})
%     axis ([-1 1 -1 1])
%     hold off
% end

%% Part III - Constructing the Julia Set
clearvars -except cval

c = -0.123 - 0.745*(1i); % need to specify value of c
M =[];  

for i=1:100                
    y = -1 + (i-1)*.01;        
    for j=1:200             
        x = -1 + (j-1)*.01;     
        z = x + 1i*y;             
        count = 0;
        while count < 100
            count = count+1;
            z = z - c;
            a = real(z);
            b = imag(z);
            r = sqrt(a^2 + b^2);
            theta = atan2(b,a);
            z = sqrt(r)*cos(theta/2)+sqrt(r)*(1i)*sin(theta/2);
            m = randi(50,1,1);
            if mod(m,2) == 1
                z = z;
            else
                z = -z;
            end;
        end;
    M = [M,z]; 
    end
end

% plot(M,'o')
% pbaspect([1 1 1]);
% axis([-1.5 1.5 -1.5 1.5])      

%% Part IV - Computing the Fractal Dimension
clearvars -except cval

phi = @(z,cval) z^2 + cval;
range = 2; 
x = linspace(-range, range, 100);
y = x;
N= 0;

for i = 1:length(cval)
    M{i} = ones(length(x), length(y));
    for j = 1:length(x)
        for k = 1:length(y)
            clear z;
            z = x(j) + 1i*y(k);
            for l = 1:100
                z(l+1) = phi( z(l), cval(i) );
                if abs( z(l+1) ) > 2
                   M{i}(j,k) = 2;
                   break;
                end
            end
        end
    end
    rval = 2*range/100; 
%     fprintf('The fractal dimension for the complex value c = %.3f + %.3f ',real(cval(i)),imag(cval(i)))
    n = length(M{i});
        for m = 1:n
            for l = 1:n
                if M{i}(m,l) == 1
                    N = N+1;
                end
            end
        end
        d = log(N)/log(1/rval);
%         fprintf('is: %5.4f\n', d); 
end

%% Part V - Connectivity of the Julia Set
clearvars -except cval 

c = -1.25; % need to specify c for each individual value
phi = @(z) z^2 - 1.25; 
z = 0; 

for k = 1:101
    z = phi(z);
    if abs(z) > 100
%         fprintf('We assume that divergence occurs for |z| > 100.\n')
%         fprintf('The orbit diverged after %.0f iterations.\n', k)
%         fprintf('We can conclude that the set is not connected for c = %.3f.\n',c)
        break;
    end
end

if abs(z) < 100
%     fprintf('We assume that divergence occurs for |z| > 100.\n')
%     fprintf('After %.0f iterations, the set failed to diverge.\n', 101)
%     fprintf('We can assume that the Julia set is connected for c = %.3f.\n',c)
end

%% Part VI: Coloring Divergent Orbits
clearvars -except cval

phi = @(z,c) z^2+ cval;
xvals = [2 2 2 2]; 
yvals = xvals; 

phi = @(z,c) z^2+ c;
for i = 1:length(cval)
    a = linspace(-xvals(i), xvals(i), 400);
    b = linspace(-yvals(i), yvals(i), 400);
    for j = 1:length(cval)
        M{j} = ones(length(a),length(b));
        for k = 1:length(a)
            for l = 1:length(b)
                clear z;
                z = a(k) + (1i)*b(l);
                for m = 1:100
                    z(m+1) = phi(z(m),cval(j));
                    if abs(z(m+1))> 100
                            M{j}(k,l) = m;
                        break;
                    end
                end
            end
        end
    end
end

% for i = 1:length(cval)
%     M{i} = M{i}';
%     figure(); hold on
%     image( [-xvals(i) xvals(i)], [-yvals(i) yvals(i)], M{i})
%     axis equal;
%     axis ([-2 2 -2 2])
%     colormap(hot(max(max(M{i})))); hold off
% end

%% Part VII: Newton's Method in the Complex Plane

phi = @(z,n) (z^n - 1)/(n*z^(n-1));
x = linspace(-5,5,600); y = x;
maxval = 5; 
M = cell(maxval-1,1);

for i = 2:maxval
    M{i-1} = 100*ones(length(x),length(y)); 
    phi = @(z) (z^i - 1)/(i*z^(i-1)); 
    for j = 1:length(x)
        for k = 1:length(y)
            z = x(j) + 1i*y(k);
            for l = 1:100
                if abs(z^i-1) > 0.001 
                    z = z - phi(z); 
                else
                    if l < 3 
%                         fprintf('If n = %i then the root of the polynomial is near (%.4f, %.4f).\n',i,x(j),y(k));
                    end
                    M{i-1}(j,k) = l; 
                    break;
                end
            end
        end
    end
end

% for i = 1:maxval-1
%     figure(); 
%     hold on
%     image( [min(x) max(x)], [min(y) max(y)], M{i})
%     axis equal
%     axis([-5 5 -5 5])
%     colormap(hot(max(max(M{i}))));
%     hold off
% end
%% Part VIII: Mandelbrot Set
clearvars -except cval

phi = @(z,cval) z^2 + cval; 
x = linspace(-2,2,500); y = x;
M = ones(length(x),length(y));

for i = 1:length(x)
    for j = 1:length(y)
        clear z;
        z = 0; 
        f = x(i) + y(j)*1i;
        for k = 1:100
            z(k+1) = phi(z(k),f);
            if abs(z(k+1)) > 100
                M(i,j) = k;
                break;
            end
        end
    end
end

% figure(); hold on
% colormap hot
% image( [-2 2], [-2 2], M')
% axis equal
% axis([ -2 1 -1.5 1.5])
