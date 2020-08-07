clc; clear; close;

n = 7; % DOF

% symbolic variables
q  = sym('q',  [n 1], 'real'); % generalized coordinates (joint angles)
qd = sym('qd', [n 1], 'real'); % "q dot" - the first derivative of the q's in time (joint velocities)
qdd = sym('qdd', [n 1], 'real');
d  = sym('d',  [n 1], 'real'); % link offsets
m  = sym('m',  [n 1], 'real'); % mass of each link
syms a1 g

% Generate inertia tensor for each link (these are wrt the inertial frame)
% Generate basic Inertia tensor
I = cell(n,1); % Inertia tensor cell
suffix = ["Ixx" "Ixy" "Ixz";...
          "Iyx" "Iyy" "Iyz";...
          "Izx" "Izy" "Izz"];
% Add index of each link
for i=1:n
    I{i,1} = sym(suffix+i);
end

% the center of mass of each link measured relative to the link fixed frame
% like Ti and Jw, c is an nx1 cell array where each element is a symoblic vector/matrix
% for example: c{3} = [c3x c3y c3z]' is the center of mass of link 3 measured relative to frame 3
c = arrayfun(@(x) [sym(['c' num2str(x) 'x'], 'real'), sym(['c' num2str(x) 'y'], 'real'), ...
    sym(['c' num2str(x) 'z'], 'real')]', 1:n, 'UniformOutput', 0)';

% cell array of your homogeneous transformations; each Ti{i} is a 4x4 symbolic transform matrix
% provide your answer in terms of the given symbolic variables
% NOTE: for symbolic arrays: q(1) = q1, q(2) = q2, etc.
Ti = cell(n,1); 

% Calculate Homogeneous transformation matrix
% Create DH parameter table
th   = [q(1,1) q(2,1) q(3,1) q(4,1) q(5,1) q(6,1) q(7,1)]; % alpha 
alp  = [-pi/2 -pi/2 -pi/2 -pi/2 -pi/2 -pi/2 0]; % theta
r    = [a1 0 0 0 0 0 0]; % r
d_    = [d(1,1) d(2,1) d(3,1) d(4,1) d(5,1) d(6,1) d(7,1)]; % d

i = 1;
Ti{i} = [cos(th(1,i)) -sin(th(1,i))*cos(alp(1,i)) sin(th(1,i))*sin(alp(1,i))  r(1,i)*cos(th(1,i));...
         sin(th(1,i)) cos(th(1,i))*cos(alp(1,i))  -cos(th(1,i))*sin(alp(1,i)) r(1,i)*sin(th(1,i));...
         0            sin(alp(1,i))               cos(alp(1,i))               d_(1,i);...
         0            0                           0                           1];

for i = 2:n
    Ti{i} = Ti{i-1} * [cos(th(1,i)) -sin(th(1,i))*cos(alp(1,i)) sin(th(1,i))*sin(alp(1,i))  r(1,i)*cos(th(1,i));...
                       sin(th(1,i)) cos(th(1,i))*cos(alp(1,i))  -cos(th(1,i))*sin(alp(1,i)) r(1,i)*sin(th(1,i));...
                       0            sin(alp(1,i))               cos(alp(1,i))               d_(1,i);...
                       0            0                           0                           1];
end

% Calculate angular velocity Jacobian 
Jw = arrayfun(@(x) sym(zeros(3,n)), 1:n, 'UniformOutput', 0)';

for i = 1:n
    Jw{i} = [[0;0;1]                      (i>1)*Ti{1}(1:3,1:3)*[0;0;1]...
             (i>2)*Ti{2}(1:3,1:3)*[0;0;1] (i>3)*Ti{3}(1:3,1:3)*[0;0;1]...
             (i>4)*Ti{4}(1:3,1:3)*[0;0;1] (i>5)*Ti{5}(1:3,1:3)*[0;0;1]...
             (i>6)*Ti{6}(1:3,1:3)*[0;0;1]];
end

% Calculate linear velocity Jacobian
Jv = cell(n,1);
R00 = eye(3);
d00 = [0;0;0];

for i = 1:n
    Jv{i} = [
        cross(R00*[0;0;1],Ti{i}(1:3,4)-d00+Ti{i}(1:3,1:3)*[c{i}(1,1);c{i}(2,1);c{i}(3,1)])... 
        (i>1)*cross(Ti{1}(1:3,1:3)*[0;0;1],Ti{i}(1:3,4)-Ti{1}(1:3,4)+Ti{i}(1:3,1:3)*[c{i}(1,1);c{i}(2,1);c{i}(3,1)])...
        (i>2)*cross(Ti{2}(1:3,1:3)*[0;0;1],Ti{i}(1:3,4)-Ti{2}(1:3,4)+Ti{i}(1:3,1:3)*[c{i}(1,1);c{i}(2,1);c{i}(3,1)])...
        (i>3)*cross(Ti{3}(1:3,1:3)*[0;0;1],Ti{i}(1:3,4)-Ti{3}(1:3,4)+Ti{i}(1:3,1:3)*[c{i}(1,1);c{i}(2,1);c{i}(3,1)])...
        (i>4)*cross(Ti{4}(1:3,1:3)*[0;0;1],Ti{i}(1:3,4)-Ti{4}(1:3,4)+Ti{i}(1:3,1:3)*[c{i}(1,1);c{i}(2,1);c{i}(3,1)])...
        (i>5)*cross(Ti{5}(1:3,1:3)*[0;0;1],Ti{i}(1:3,4)-Ti{5}(1:3,4)+Ti{i}(1:3,1:3)*[c{i}(1,1);c{i}(2,1);c{i}(3,1)])...
        (i>6)*cross(Ti{6}(1:3,1:3)*[0;0;1],Ti{i}(1:3,4)-Ti{6}(1:3,4)+Ti{i}(1:3,1:3)*[c{i}(1,1);c{i}(2,1);c{i}(3,1)])];
end

% Calculate potential energy
g_ = [0 0 g]; % gravity vector
PE = 0;

for i=1:n
    dist = Ti{i}*[c{i};1];
    dist = dist(1:3,1);
    PE = PE + m(i,1)*g_*dist;
end

% Calculate kinetic energy
D = 0;
for i = 1:n
    D = D + m(i,1)*Jv{i}'*Jv{i} + Jw{i}'*I{i}*Jw{i};
end

KE = 0.5 * qd' * D * qd;

% Calculate coriolis and centrifugal forces
C = sym(zeros(n));

for j = 1:n
    for k = 1:n
        for i = 1:n
            C(j,k) = C(j,k) + 0.5*( diff(D(k,j),q(i)) + diff(D(k,i),q(j)) - diff(D(i,j),q(k)) ) * qd(i);
        end
    end
end

eom_lhs = D * qdd + C * qd + [diff(PE, q(1)); diff(PE, q(2)); diff(PE, q(3));...
                              diff(PE, q(4)); diff(PE, q(5)); diff(PE, q(6)); diff(PE, q(7))];