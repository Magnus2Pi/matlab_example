% FUNCTION Newton-Raphson

% INPUT
% Three functions K, mu, rho : D \subset R3 -> R
% Three values K_g, mu_g, rho_g 

% OUTPUT 
% { p = p(x,y,z) \elem D | K(p) = K_g  and  mu(p) = mu_g  and  rho(p) = rho_g }

% INPUT
% xc is (equidistant) discretisation of the x-axis
% yc is (equidistant) discretisation of the y-axis
% zc is (equidistant) discretisation of the z-axis
% K, mu and rho are functions of (xc,yc,zc)
% K_g, mu_g and rho_g are given (measured) values
% x,y,z can be porosity, lithology, saturation
% K, mu, rho can be bulk modulus, shear modulus, density

% OUTPUT
% intersection_points that are all the common intersection points of the
% equisurfaces K = K_g, mu = mu_g and rho = rho_g

% If we draw the equisurfaces K = K_g etc. the solution to the inverse
% problem is all the common intersections of the three surfaces.
% The problem is that the functions K etc. is only given in the
% discrete points (xc,yc,zc)

% The algorithm works as follows:
% 1) First a search through all the small cubes is made since:
% A necessary condition for the solution to be in the small cube
% (xc_i, xc_i+1) X (yc_j, yc_j+1) X (zc_k,zc_k+1) is that
% the values for K_g etc. lies in the interval [K_min , K_max],
% where K_min is the minimum K-value among the eight corners of the small cube.

% 2) For each small cube (where the solution might be) the given discrete
% functions K etc. are interpolated using three-linear shape functions.
% The interpolating functions are polynomials in x,y,z which are symbolic
% so they can be operated on by maple. (e.g. differentiated)

% 3) To find the common intersection of the three surfaces Newton-Raphson
% method generalised for 3D is used. (In 1D: x_n+1 = x_n - f(x_n)/f´(x_n))
% This is an iterative method, so we need an initial value which is set to 
% the center of the cube. 

% 4) Only solutions which are inside the small cube under consideration are
% kept as true solutions. 

function [intersection_points] = newtonR(xc,yc,zc,rho,K,mu,rho_g,K_g,mu_g)
[ny nx nz] = size(rho);
dx         = xc(2)-xc(1);
dy         = yc(2)-yc(1);
dz         = zc(2)-zc(1);

% 1)
% the "possible_solution" array contains after this all possible cubes where point of intersection might be
possible_solutions = [];
for i=1:nx-1 % x
    for j=1:ny-1 % y
        for k=1:nz-1 % z
            rho_min = min( ...
                [rho(j  ,i  ,k),rho(j+1,i,k  ),rho(j,i+1,k  ),rho(j  ,i  ,k+1) ...
                ,rho(j+1,i+1,k),rho(j+1,i,k+1),rho(j,i+1,k+1),rho(j+1,i+1,k+1)]);
            rho_max = max( ...
                [rho(j  ,i  ,k),rho(j+1,i,k  ),rho(j,i+1,k  ),rho(j  ,i  ,k+1) ...
                ,rho(j+1,i+1,k),rho(j+1,i,k+1),rho(j,i+1,k+1),rho(j+1,i+1,k+1)]);
            K_min   = min( ...
                [K(j  ,i  ,k),  K(j+1,i,k  ),  K(j,i+1,k  ),  K(j  ,i  ,k+1) ...
                ,K(j+1,i+1,k),  K(j+1,i,k+1),  K(j,i+1,k+1),  K(j+1,i+1,k+1)]);
            K_max   = max( ...
                [K(j  ,i  ,k),  K(j+1,i,k  ),  K(j,i+1,k  ),  K(j  ,i  ,k+1) ...
                ,K(j+1,i+1,k),  K(j+1,i,k+1),  K(j,i+1,k+1),  K(j+1,i+1,k+1)]);
            mu_min  = min( ...
                [mu(j  ,i  ,k), mu(j+1,i,k  ), mu(j,i+1,k  ), mu(j  ,i  ,k+1) ...
                ,mu(j+1,i+1,k), mu(j+1,i,k+1), mu(j,i+1,k+1), mu(j+1,i+1,k+1)]);
            mu_max  = max( ...
                [mu(j  ,i  ,k), mu(j+1,i,k  ), mu(j,i+1,k  ), mu(j  ,i  ,k+1) ...
                ,mu(j+1,i+1,k), mu(j+1,i,k+1), mu(j,i+1,k+1), mu(j+1,i+1,k+1)]);
            if(     rho_min <= rho_g && rho_g <= rho_max && ...
                    K_min   <= K_g   && K_g   <= K_max   && ...
                    mu_min  <= mu_g  && mu_g  <= mu_max)
                possible_solutions = [possible_solutions ; i j k ];
            end
        end
    end
end

% 2)
syms x real;
syms y real;
syms z real;

num_possible_solutions = size(possible_solutions,1);
for c = 1:num_possible_solutions
    cube_ind  = possible_solutions(c,:);
    i = cube_ind(1);
    j = cube_ind(2);
    k = cube_ind(3);
    
    V    = (xc(i+1)-xc(i))*(yc(j+1)-yc(j))*(zc(k+1)-zc(k));
    %            y  ,x  ,z
    % value at each corner
    rho_1  = rho(j  ,i  ,k  );
    rho_2  = rho(j  ,i+1,k  );
    rho_3  = rho(j+1,i+1,k  );
    rho_4  = rho(j+1,i  ,k  );
    rho_5  = rho(j  ,i  ,k+1);
    rho_6  = rho(j  ,i+1,k+1);
    rho_7  = rho(j+1,i+1,k+1);
    rho_8  = rho(j+1,i  ,k+1);
    % shape function for each corner point
    s1_rho = rho_1*(xc(i+1)-x    )*(yc(j+1)-y    )*(zc(k+1)-z    )/V;
    s2_rho = rho_2*(x      -xc(i))*(yc(j+1)-y    )*(zc(k+1)-z    )/V;
    s3_rho = rho_3*(x      -xc(i))*(y      -yc(j))*(zc(k+1)-z    )/V;
    s4_rho = rho_4*(xc(i+1)-x    )*(y      -yc(j))*(zc(k+1)-z    )/V;
    s5_rho = rho_5*(xc(i+1)-x    )*(yc(j+1)-y    )*(z      -zc(k))/V;
    s6_rho = rho_6*(x      -xc(i))*(yc(j+1)-y    )*(z      -zc(k))/V;
    s7_rho = rho_7*(x      -xc(i))*(y      -yc(j))*(z      -zc(k))/V;
    s8_rho = rho_8*(xc(i+1)-x    )*(y      -yc(j))*(z      -zc(k))/V;
    % the interpolating function now have the right values in the corners
    rho_i  = simplify(s1_rho + s2_rho + s3_rho + s4_rho + s5_rho + s6_rho + s7_rho + s8_rho);
    
    %        y  ,x  ,z
    K_1  = K(j  ,i  ,k  );
    K_2  = K(j  ,i+1,k  );
    K_3  = K(j+1,i+1,k  );
    K_4  = K(j+1,i  ,k  );
    K_5  = K(j  ,i  ,k+1);
    K_6  = K(j  ,i+1,k+1);
    K_7  = K(j+1,i+1,k+1);
    K_8  = K(j+1,i  ,k+1);
    s1_K = K_1*(xc(i+1)-x    )*(yc(j+1)-y    )*(zc(k+1)-z    )/V;
    s2_K = K_2*(x      -xc(i))*(yc(j+1)-y    )*(zc(k+1)-z    )/V;
    s3_K = K_3*(x      -xc(i))*(y      -yc(j))*(zc(k+1)-z    )/V;
    s4_K = K_4*(xc(i+1)-x    )*(y      -yc(j))*(zc(k+1)-z    )/V;
    s5_K = K_5*(xc(i+1)-x    )*(yc(j+1)-y    )*(z      -zc(k))/V;
    s6_K = K_6*(x      -xc(i))*(yc(j+1)-y    )*(z      -zc(k))/V;
    s7_K = K_7*(x      -xc(i))*(y      -yc(j))*(z      -zc(k))/V;
    s8_K = K_8*(xc(i+1)-x    )*(y      -yc(j))*(z      -zc(k))/V;
    K_i  = simplify(s1_K + s2_K + s3_K + s4_K + s5_K + s6_K + s7_K + s8_K);
    
    %          y  ,x  ,z
    mu_1  = mu(j  ,i  ,k  );
    mu_2  = mu(j  ,i+1,k  );
    mu_3  = mu(j+1,i+1,k  );
    mu_4  = mu(j+1,i  ,k  );
    mu_5  = mu(j  ,i  ,k+1);
    mu_6  = mu(j  ,i+1,k+1);
    mu_7  = mu(j+1,i+1,k+1);
    mu_8  = mu(j+1,i  ,k+1);
    s1_mu = mu_1*(xc(i+1)-x    )*(yc(j+1)-y    )*(zc(k+1)-z    )/V;
    s2_mu = mu_2*(x      -xc(i))*(yc(j+1)-y    )*(zc(k+1)-z    )/V;
    s3_mu = mu_3*(x      -xc(i))*(y      -yc(j))*(zc(k+1)-z    )/V;
    s4_mu = mu_4*(xc(i+1)-x    )*(y      -yc(j))*(zc(k+1)-z    )/V;
    s5_mu = mu_5*(xc(i+1)-x    )*(yc(j+1)-y    )*(z      -zc(k))/V;
    s6_mu = mu_6*(x      -xc(i))*(yc(j+1)-y    )*(z      -zc(k))/V;
    s7_mu = mu_7*(x      -xc(i))*(y      -yc(j))*(z      -zc(k))/V;
    s8_mu = mu_8*(xc(i+1)-x    )*(y      -yc(j))*(z      -zc(k))/V;
    mu_i  = simplify(s1_mu + s2_mu + s3_mu + s4_mu + s5_mu + s6_mu + s7_mu + s8_mu);
    
    rho_K_mu{c}   = {rho_i,K_i,mu_i};
end

% 3)
nit     = 10;       % max iterations
epsilon = 1.0e-12;  % stopp interating if distance is less that epsilon
intersection_points   = [];
for c = 1:num_possible_solutions % for each cube
    % the interpolating functions
    rho_i = rho_K_mu{c}{1};
    K_i   = rho_K_mu{c}{2};
    mu_i  = rho_K_mu{c}{3};
    
    cube_ind  = possible_solutions(c,:);
    i = cube_ind(1);
    j = cube_ind(2);
    k = cube_ind(3);
    cord_corner = [xc(i) yc(j) zc(k)];
    cord_center = cord_corner + 1/2*[dx dy dz];  
    % initial guess
    x = cord_center(1);
    y = cord_center(2);
    z = cord_center(3);
    
    ni  =  0;
    delta = 2*epsilon;
    % Newton-Raphson in 3D
    while delta > epsilon && ni < nit
        F(1,1) = eval( rho_i - rho_g);
        F(2,1) = eval( K_i   - K_g  );
        F(3,1) = eval( mu_i  - mu_g );
        
        A(1,1) = eval( diff(rho_i,'x') );
        A(1,2) = eval( diff(rho_i,'y') );
        A(1,3) = eval( diff(rho_i,'z') );
        
        A(2,1) = eval( diff(K_i  ,'x') );
        A(2,2) = eval( diff(K_i  ,'y') );
        A(2,3) = eval( diff(K_i  ,'z') );
        
        A(3,1) = eval( diff(mu_i ,'x') );
        A(3,2) = eval( diff(mu_i ,'y') );
        A(3,3) = eval( diff(mu_i ,'z') );
        
        b      = [x y z]' - A\F;
        
        delta = norm(b - [x ; y ; z]);
        
        x      = b(1);
        y      = b(2);
        z      = b(3);
        
        ni     = ni + 1;
    end % while
    
    % 4)
    % only keep solutions that are inside the small cube
    if (cord_corner(1) <= x && x <= cord_corner(1)+dx && ...
            cord_corner(2) <= y && y <= cord_corner(2)+dy && ...
            cord_corner(3) <= z && z <= cord_corner(3)+dz)
        intersection_points = [intersection_points ; [x y z c]];
    end % if    
end % for

end % function


