% Define the function "plotsolution_HW6" with inputs A and b,
% and provide a comment that describes what the function does
% and how it is used.

function [answer] = plotsolution_HW6(A, b)

%PLOTSOLUTION_HW6 receives A and B and plots the particular solution and
%spanning vectors in the solution


% Assign the number of rows of "A" to a variable "m", 
% and the number of columns to the variable "n".
% Check that n is 2 or 3. Otherwise, give an error message.

m = size(A, 1);
n = size(A, 2);

if n ~= 2 && n ~= 3
    error("The matrix must have 2 or 3 columns")
end


% Include this for numerical stability.
if m == n
    A = [A;zeros(1,n)];
    b = [b;0];
end

% Check for inconsistency. (Hint: use rref.) 
% If the system is consistent, then find a particular solution p.

[~,piv] = rref([A b]);

if piv(end)==(n+1)
    error('The system of equations is inconsistent!');    
else
    p = A\b;
end 

% Find vectors that span the set of homogeneous solutions.
V = null(A,'r');
num_vectors = size(V,2); 

% This creates the figure and orients the plot.
% You can modify the limits to change the point of view.
figure; hold on;
xlim([p(1)-2 p(1)+2]);ylim([p(2)-2 p(2)+2]);
if n == 3
    zlim([p(3)-2 p(3)+2])
end

% This plots the particular solution p as an arrow
if n == 3
    view([-37.5,30]);
    scatter3(0,0,0,'go','LineWidth',3);
    quiver3(0,0,0,p(1),p(2),p(3),1.0,'LineWidth',3,'MaxHeadSize',.5,'Color',[0 0 0]);
else
    plot(0,0,'go','LineWidth',3);
    quiver(0,0,p(1),p(2),1.0,'LineWidth',3,'MaxHeadSize',.5,'Color',[0 0 0]);
end
axis square
xlabel('x1');ylabel('x2');zlabel('x3')



% This plots the basis vectors as arrows.
for ii = (1:num_vectors)
    if n == 3
        quiver3(0,0,0,V(1,ii),V(2,ii),V(3,ii),1.0,'LineWidth',3,'MaxHeadSize',.5,'Color',[1 0 0]);
    else
        quiver(0,0,V(1,ii),V(2,ii),1.0,'LineWidth',3,'MaxHeadSize',.5,'Color',[1 0 0]);
    end
end
leg = legend('Origin','Particular solution','Spanning vectors','AutoUpdate','off');

num_points = 200;

% Initialize a matrix "points" of zeros with "n" rows and "num_points" columns.

points = zeros(n, num_points);

% Compute and plot the points one by one.
% Create a loop so that the number of points plotted is num_points. 
% Note that a point is discarded if it falls outside the given region 
% specified in the following "if" statement.

ii = 1;
while ii<=num_points
    
    % Inside the loop:
    % Generate a solution p + V*r where r is a random vector of combining
    % coefficients with elements between -50 and 50. 
    % Use the function rand, which generates random numbers between 0 and 1. 
    % Note that (rand - 0.5) shifts the range to -0.5 to 0.5.
    
    % First assign the vector r.
   r = -50 +(100).*rand(num_vectors,1);
    
    % This keeps the points to be plotted within a desired range.
   if abs(V*r) < 2*ones(n,1)
        % In the ii'th column of "points", store the linear combination of the
        % columns of "V" weighted by the values in "r", plus your particular solution "p".
        points(:,ii)=(V*r)+p;    
        
        % This adds the new point to the plot.
        if n == 3
            scatter3(points(1,ii),points(2,ii),points(3,ii),'b*');
        else
            scatter(points(1,ii),points(2,ii),'b*');
        end
    
        % This pauses the program so that you can see each point as it is added.
        pause(.01)
        ii=ii+1;
   end
        %Complete the loop.
end


% Needed for the legend.
h.DisplayName = 'b';


% Test Cases:
%   3D:
%   A = [1 3 2;0 0 0;-1 -2 -3]'; b = [1;2;3]; plotsolution_HW6(A,b);
%   A = [-2 2 4;-2 2 5;1 -1 -1]; b = [-6;-7;2]; plotsolution_HW6(A,b);
%   r = rand(3,1); A = [r 2*r -r]; plotsolution_HW6(A,r);
%   A = rand(3,3); b = rand(3,1); plotsolution_HW6(A,b);
%   r = rand(3,2); A = [r 2*r(:,1)]; b = rand(3,1); plotsolution_HW6(A,b);
%   r = rand(3,2); A = [r 2*r(:,1)]; b = A(:,1) - A(:,2); plotsolution_HW6(A,b);
%   c = [1 3 2]; A = [c;-2*c;3*c;6*c]; b = A(:,1)-A(:,2)-A(:,3); plotsolution_HW6(A,b);
%
%   2D:
%   A = [1 -2;-2 4]; b = [1;2]; plotsolution_HW6(A,b);
%   A = [1 -2;-2 4]; b = [3;-6]; plotsolution_HW6(A,b);
%   A = rand(2,2); b = rand(2,1); plotsolution_HW6(A,b);
%   A = rand(1,2); b = rand(1,1); plotsolution_HW6(A,b);
%   A = rand(5,2); b = rand(5,1); plotsolution_HW6(A,b);
%   r = rand(5,1); A = [r 2*r]; b = 3*r; plotsolution_HW6(A,b);
