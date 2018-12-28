function [ outT, outX, outVal, outGr, evalNumbers] = MoreThuente(functionName, params)

%   ------------------      *******************        ------------------
%   *                                                                   *
%   *               *************************************               *
%   *               *                                   *               *
%   *               *      Moré Thuente line search     *               *
%   *               *                                   *               *
%   *               *************************************               *
%   *                                                                   *
%   ------------------      *******************        ------------------

%   Moré Thuente line search is a line search procedure for computing 
%   step-size parameter t such that it satisfies so called strong 
%   Wolfe condition. It is an iterative method which is proven
%   to be very effective. The authors claim that the algorithm 
%   terminates within a small number of iterations. It is one of the 
%   most popular method in line search cathegory. This method is 
%   originally proposed by J.J. Moré and D.J. Thuente.

%   J.J. Moré, D.J. Thuente,
%   Line search algorithms with guaranteed sufficient decrease,
%   ACM Trans.Math. Softw., 20 (3) (1994) 286-307.

%   ------------------      *******************        ------------------

    evalNumbers = EvaluationNumbers(0,0,0);
    
    %----------------------------------------------------------------------
    % parameters that are set through the gui application 
    x_0 = params.startingPoint;
    g_at_x_0 = params.grad;
    d = params.dir;
    f_at_x_0 = params.vals(size(params.vals,2));
    mu = params.rho; 
    alpha_t = params.tInitStart;
    eta = params.sigma;
    workPrec = 1e-16;
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % all the vectors are of N x 1 dimension inside this method 
    if(size(x_0,2) > 1)
        x_0 = x_0';
    end
    if(size(g_at_x_0,2) > 1)
        g_at_x_0 = g_at_x_0';
    end
    if(size(d,2) > 1)
        d = d';
    end
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % some local constants 
    p50 = .5;
    p66 = .66;
    delta_max = 4;
    alpha_max0 = 100;
    alpha_min0 = 0; 
    max_iter = 50;
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % parameters needed for debugging 
    n_iter = 0;
    alpha_l = 0.0;
    alpha_u = 0.0;
    isBracketed = false;
    stage1 = true;
    width = alpha_max0 - alpha_min0;
    width1 = 2*width;
        
    f = functionName;    
    phi_prim_0 = g_at_x_0'*d;
    
    if alpha_l==alpha_u
       if alpha_l==0
           phi_psi_alpha_u = f_at_x_0;
           phi_psi_alpha_l = f_at_x_0;
           g_u=g_at_x_0;
           g_l=g_at_x_0;
           
       else
           [phi_psi_alpha_u,g_u,~] = feval(f,x_0 + alpha_u*d,[1 1 0]);
           phi_psi_alpha_l = phi_psi_alpha_u;
           g_l=g_u;
           evalNumbers.incrementBy([1 1 0]);
       end
    else
        [phi_psi_alpha_u,g_u,~] = feval(f,x_0 + alpha_u*d,[1 1 0]);
        [phi_psi_alpha_l,g_l,~] = feval(f,x_0 + alpha_l*d,[1 1 0]);
        evalNumbers.incrementBy([2 2 0]);
    end
        
    [phi_psi_alpha_t,g_t,~] = feval(f,x_0 + alpha_t*d,[1 1 0]);
    evalNumbers.incrementBy([1 1 0]);
    
    while true  
       n_iter = n_iter + 1;   
       if isBracketed
           alpha_min = min(alpha_l,alpha_u);
           alpha_max = max(alpha_l,alpha_u);
       else
           alpha_min = alpha_l;
           alpha_max = alpha_t + delta_max*(alpha_t-alpha_l);
       end
       
       if alpha_t < alpha_min0
          alpha_t = alpha_min0;
          [phi_psi_alpha_min0,g_min0,~] = feval(f,x_0 + alpha_min0*d,[1 1 0]);
          g_t=g_min0;
          phi_psi_alpha_t = phi_psi_alpha_min0;
          evalNumbers.incrementBy([1 1 0]);
       end
       if alpha_t > alpha_max0
          alpha_t = alpha_max0;
          [phi_psi_alpha_max0,g_max0,~] = feval(f,x_0 + alpha_max0*d,[1 1 0]);
          g_t=g_max0;
          phi_psi_alpha_t = phi_psi_alpha_max0;
          evalNumbers.incrementBy([1 1 0]);
       end
       
       if (isBracketed && abs(alpha_l-alpha_u) < workPrec*alpha_max)
           
           outT = alpha_l;
           outX = (x_0+alpha_l*d)';
           outVal = phi_psi_alpha_l;
           outGr = g_l;
           return;
       end
       
       if (alpha_t > 0 && phi_psi_alpha_t <= f_at_x_0 + mu*alpha_t*phi_prim_0 && ...
           abs(g_t'*d) <= eta*abs(phi_prim_0))       
           
           outT = alpha_t;
           outX = (x_0+alpha_t*d)';
           outVal = phi_psi_alpha_t;
           outGr = g_t;
           return;
       end
           
       if (stage1 && phi_psi_alpha_t - f_at_x_0 - mu*alpha_t*phi_prim_0 <= 0 && g_t'*d - mu*phi_prim_0 > 0)
          stage1 = false;
       end
      
       if(stage1 && phi_psi_alpha_t <= phi_psi_alpha_l &&  phi_psi_alpha_t > f_at_x_0 + mu*alpha_t*phi_prim_0)
           c1 = f_at_x_0 + mu*alpha_t*phi_prim_0;
           c2 = mu*phi_prim_0;           
           [alpha_t_new,isBracketed,refinePoint]=nextSafeguard2(alpha_u,alpha_l,alpha_t,phi_psi_alpha_u-c1,phi_psi_alpha_l-c1,phi_psi_alpha_t-c1,g_u,g_l,g_t,isBracketed,workPrec,alpha_max,alpha_min,d,c2);
           [alpha_u,alpha_l,phi_psi_alpha_u,phi_psi_alpha_l,g_u,g_l] = generateNextInterval2(alpha_u,alpha_l,alpha_t,phi_psi_alpha_u-c1,phi_psi_alpha_l-c1,phi_psi_alpha_t-c1,g_u,g_l,g_t,d);
       else
           [alpha_t_new,isBracketed,refinePoint]=nextSafeguard2(alpha_u,alpha_l,alpha_t,phi_psi_alpha_u,phi_psi_alpha_l,phi_psi_alpha_t,g_u,g_l,g_t,isBracketed,workPrec,alpha_max,alpha_min,d,0);
           [alpha_u,alpha_l,phi_psi_alpha_u,phi_psi_alpha_l,g_u,g_l] = generateNextInterval2(alpha_u,alpha_l,alpha_t,phi_psi_alpha_u,phi_psi_alpha_l,phi_psi_alpha_t,g_u,g_l,g_t,d);
       end
       
       alpha_t = alpha_t_new;
       
       if alpha_t > alpha_max
           alpha_t = alpha_max;
       end
       if alpha_t < alpha_min
           alpha_t = alpha_min; 
       end
        
       if (isBracketed && refinePoint)
           if (alpha_u > alpha_l) 
             alpha_t = min(alpha_l+p66*(alpha_u-alpha_l),alpha_t);
          else
             alpha_t = max(alpha_l+p66*(alpha_u-alpha_l),alpha_t);
          end
       end      
       
       if isBracketed 
           if (abs(alpha_l-alpha_u) >= p66*width1)
               alpha_t = p50*(alpha_l + alpha_u);
           end
           width1 = width;
           width = abs(alpha_l-alpha_u);
       end
           
       if (n_iter == max_iter)
           
           outT = alpha_l;
           outX = (x_0+alpha_l*d)';           
           outVal = phi_psi_alpha_l;
           outGr = g_l;
           return;
       end
       
       [phi_psi_alpha_t,g_t,~] = feval(f,x_0 + alpha_t*d,[1 1 0]);
       evalNumbers.incrementBy([1 1 0]);
    end

    function [alpha_t_new,isBracketed,refinePoint] = nextSafeguard2(alpha_u,alpha_l,alpha_t,phi_psi_alpha_u,phi_psi_alpha_l,phi_psi_alpha_t,g_u,g_l,g_t,isBracketed,workPrec,alpha_max,alpha_min,d,c2)
        phi_psi_prim_alpha_u_loc = g_u'*d-c2;
        phi_psi_prim_alpha_l_loc = g_l'*d-c2;
        phi_psi_prim_alpha_t_loc = g_t'*d-c2;
        
        if(phi_psi_alpha_t > phi_psi_alpha_l)
            
            if(alpha_t == alpha_l)                
                theta = 3.*(phi_psi_alpha_l - phi_psi_alpha_t)/(alpha_t - alpha_l + workPrec) + phi_psi_prim_alpha_l_loc + phi_psi_prim_alpha_t_loc;
            else                
                theta = 3.*(phi_psi_alpha_l - phi_psi_alpha_t)/(alpha_t - alpha_l) + phi_psi_prim_alpha_l_loc + phi_psi_prim_alpha_t_loc;
            end
            
            s = norm([theta,phi_psi_prim_alpha_l_loc,phi_psi_prim_alpha_t_loc],inf);
            gamma = s*sqrt((theta/s).^2 - (phi_psi_prim_alpha_l_loc/s).*(phi_psi_prim_alpha_t_loc/s));
            
            if (alpha_t < alpha_l) 
                gamma = -gamma;
            end
            
            p = (gamma - phi_psi_prim_alpha_l_loc) + theta;
            q = ((gamma - phi_psi_prim_alpha_l_loc) + gamma) + phi_psi_prim_alpha_t_loc;
            r = p/q;
            alpha_c = alpha_l + r*(alpha_t - alpha_l);
            alpha_q = alpha_l + ((phi_psi_prim_alpha_l_loc/((phi_psi_alpha_l - phi_psi_alpha_t)/(alpha_t-alpha_l)+phi_psi_prim_alpha_l_loc))/2).*(alpha_t - alpha_l);
            
            if abs(alpha_c-alpha_l) < abs(alpha_q-alpha_l)
                alpha_t_new = alpha_c;
            else
                alpha_t_new = 0.5*(alpha_q+alpha_c);
            end
            
            isBracketed = true;
            refinePoint = true;                       
        else
            
            if(alpha_t==alpha_l)
                theta = 3*(phi_psi_alpha_l - phi_psi_alpha_t)/(alpha_t - alpha_l + workPrec) + phi_psi_prim_alpha_l_loc + phi_psi_prim_alpha_t_loc;
            else               
                theta = 3*(phi_psi_alpha_l - phi_psi_alpha_t)/(alpha_t - alpha_l) + phi_psi_prim_alpha_l_loc + phi_psi_prim_alpha_t_loc; 
            end
            
            s = norm([theta,phi_psi_prim_alpha_l_loc,phi_psi_prim_alpha_t_loc],inf);
            gamma = s*sqrt((theta/s).^2 - (phi_psi_prim_alpha_l_loc/s).*(phi_psi_prim_alpha_t_loc/s));
            
            if (alpha_t > alpha_l) 
                gamma = -gamma;
            end

            p = (gamma - phi_psi_prim_alpha_t_loc) + theta;
            q = ((gamma - phi_psi_prim_alpha_t_loc) + gamma) + phi_psi_prim_alpha_l_loc;
            r = p/q;
            alpha_c = alpha_t + r*(alpha_l - alpha_t);
            alpha_s = alpha_t + (phi_psi_prim_alpha_t_loc/(phi_psi_prim_alpha_t_loc-phi_psi_prim_alpha_l_loc))*(alpha_l - alpha_t);   
            
            if phi_psi_prim_alpha_t_loc*phi_psi_prim_alpha_l_loc < 0
                
                if abs(alpha_c-alpha_t) > abs(alpha_s-alpha_t)
                   alpha_t_new = alpha_c;
                else
                   alpha_t_new = alpha_s;
                end    
                
                isBracketed = true;
                refinePoint = false;           
            else
                if abs(phi_psi_prim_alpha_t_loc) < abs(phi_psi_prim_alpha_l_loc)
                    
                    if alpha_t == alpha_l
                        theta = 3*(phi_psi_alpha_l - phi_psi_alpha_t)/(alpha_t - alpha_l+workPrec) + phi_psi_prim_alpha_l_loc + phi_psi_prim_alpha_t_loc;
                    else                       
                        theta = 3*(phi_psi_alpha_l - phi_psi_alpha_t)/(alpha_t - alpha_l) + phi_psi_prim_alpha_l_loc + phi_psi_prim_alpha_t_loc; 
                    end
                    
                    s = norm([theta,phi_psi_prim_alpha_l_loc,phi_psi_prim_alpha_t_loc],inf);
                    gamma = s*sqrt(max(0.,(theta/s)^2 - (phi_psi_prim_alpha_l_loc/s)*(phi_psi_prim_alpha_t_loc/s)));
                    
                    if (alpha_t > alpha_l) 
                        gamma = -gamma;
                    end
                    
                    p = (gamma - phi_psi_prim_alpha_t_loc) + theta;
                    q = (gamma + (phi_psi_prim_alpha_l_loc - phi_psi_prim_alpha_t_loc)) + gamma;
                    r = p/q;
                    
                    if(r < 0.0 && gamma ~= 0.0)
                       alpha_c = alpha_t + r*(alpha_l-alpha_t); 
                    elseif (alpha_t > alpha_l)
                       alpha_c = alpha_max;
                    else
                       alpha_c = alpha_min;
                    end
                    
                    if (isBracketed)
                       if (abs(alpha_t-alpha_c) < abs(alpha_t-alpha_s)) 
                          alpha_t_new = alpha_c;
                       else
                          alpha_t_new = alpha_s;
                       end
                    else
                       if (abs(alpha_t-alpha_c) > abs(alpha_t-alpha_s)) 
                          alpha_t_new = alpha_c;
                       else
                          alpha_t_new = alpha_s;
                       end 
                    end 
                    refinePoint = true;         
                else
                    if isBracketed
                        
                        if alpha_u==alpha_t
                            theta = 3*(phi_psi_alpha_t - phi_psi_alpha_u)/(alpha_u - alpha_t+workPrec) + phi_psi_prim_alpha_u_loc + phi_psi_prim_alpha_t_loc;
                        else
                            theta = 3*(phi_psi_alpha_t - phi_psi_alpha_u)/(alpha_u - alpha_t) + phi_psi_prim_alpha_u_loc + phi_psi_prim_alpha_t_loc;
                        end
                        
                        s = norm([theta,phi_psi_prim_alpha_u_loc,phi_psi_prim_alpha_t_loc],inf);
                        gamma = s*sqrt((theta/s)^2 - (phi_psi_prim_alpha_u_loc/s)*(phi_psi_prim_alpha_t_loc/s));
                        
                        if (alpha_t > alpha_u) 
                            gamma = -gamma;
                        end
                        
                        p = (gamma - phi_psi_prim_alpha_t_loc) + theta;
                        q = ((gamma - phi_psi_prim_alpha_t_loc) + gamma) + phi_psi_prim_alpha_u_loc;
                        r = p/q;
                        alpha_t_new = alpha_t + r*(alpha_u - alpha_t);
                        
                     elseif (alpha_t > alpha_l)
                        alpha_t_new = alpha_max;
                     else
                        alpha_t_new = alpha_min;
                    end
                    
                    refinePoint = false;
                end
            end
        end
    end

    function [alpha_u_new,alpha_l_new,phi_psi_alpha_u_new,phi_psi_alpha_l_new,g_u_new,g_l_new] = generateNextInterval2(alpha_u,alpha_l,alpha_t,phi_psi_alpha_u,phi_psi_alpha_l,phi_psi_alpha_t,g_u,g_l,g_t,d)
       
       alpha_u_new=alpha_u;
       alpha_l_new=alpha_l;
       phi_psi_alpha_u_new=phi_psi_alpha_u;
       phi_psi_alpha_l_new=phi_psi_alpha_l;
       g_u_new = g_u;
       g_l_new = g_l;
       
       if (phi_psi_alpha_t > phi_psi_alpha_l)
           alpha_u_new = alpha_t;
           phi_psi_alpha_u_new = phi_psi_alpha_t;
           g_u_new = g_t;
       else
           v = g_t'*d*(alpha_t-alpha_l);
           if (v < 0)
               alpha_l_new = alpha_t;
               phi_psi_alpha_l_new = phi_psi_alpha_t;
               g_l_new = g_t;
           else
               alpha_u_new = alpha_l;
               alpha_l_new = alpha_t;
               phi_psi_alpha_u_new = phi_psi_alpha_l;
               g_u_new = g_l;
               phi_psi_alpha_l_new = phi_psi_alpha_t;
               g_l_new = g_t;
           end
       end
    end
        
end
