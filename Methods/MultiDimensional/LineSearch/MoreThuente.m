%function [ outT, outX, outVal, outGr, evalNumbers, n_iter, termCrit ] = MoreThuente(functionName, params) 
function [ outT, outX, outVal, outGr, evalNumbers] = MoreThuente(functionName, params) 

    evalNumbers = EvaluationNumbers(0,0,0);
    
    %----------------------------------------------------------------------
    % parametri koji su podeseni u gui-u aplikacije
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
    % svi vektori su dimenzija N x 1 dok radim s njima u ovoj funkciji
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
    % konstanti parametri algoritma   
    p50 = .5;
    p66 = .66;
    delta_max = 4;
    alpha_max0 = 100;
    alpha_min0 = 0;  
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % parametri za pracenje algoritma
    n_iter = 0;
    alpha_l = 0.0;
    alpha_u = 0.0;
    isBracketed = false;
    refinePoint = false;
    stage1 = true;
    width = alpha_max0 - alpha_min0;
    width1 = 2*width;
    
        %------------------------------------------------------------------
        % parametar kojim se odredjuje razlog zavrsetka algoritma
        % termCrit = 0;
        % termCrit = 1 nadjena je tacka koja zadovoljava uslove
        % termCrit = 2 broj iteracija je dostigo maxIter
        % termCrit = 3 interval je previse suzen
        %------------------------------------------------------------------
        
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % moji parametri
    max_iter = 50;
    %----------------------------------------------------------------------

    f = functionName;    
    phi_prim_0 = g_at_x_0'*d;
    
    %if abs(alpha_l - alpha_u) < workPrec
    if alpha_l==alpha_u
       %if abs(alpha_l)<workPrec
       if alpha_l==0
           phi_psi_alpha_u = f_at_x_0;
           phi_psi_alpha_l = f_at_x_0;
           g_1 = g_at_x_0;
           phi_psi_prim_alpha_u = g_at_x_0'*d;
           phi_psi_prim_alpha_l = g_at_x_0'*d;
       else
           [phi_psi_alpha_u,g_u,~] = feval(f,x_0 + alpha_u*d,[1 1 0]);
           phi_psi_prim_alpha_u=g_u'*d;
           phi_psi_alpha_l = phi_psi_alpha_u;
           g_1 = g_u;
           phi_psi_prim_alpha_l = phi_psi_prim_alpha_u;
           evalNumbers.incrementBy([1 1 0]);
       end
    else
        [phi_psi_alpha_u,g_u,~] = feval(f,x_0 + alpha_u*d,[1 1 0]);
        phi_psi_prim_alpha_u=g_u'*d;
        [phi_psi_alpha_l,g_l,~] = feval(f,x_0 + alpha_l*d,[1 1 0]);
        phi_psi_prim_alpha_l=g_l'*d;
        evalNumbers.incrementBy([2 2 0]);
    end
    
    
    [phi_psi_alpha_t,g_t,~] = feval(f,x_0 + alpha_t*d,[1 1 0]);
    phi_psi_prim_alpha_t=g_t'*d;
    
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
          phi_psi_prim_alpha_min0=g_min0'*d;
          phi_psi_prim_alpha_t = phi_psi_prim_alpha_min0;
          phi_psi_alpha_t = phi_psi_alpha_min0;
          g_t = g_min0;
          evalNumbers.incrementBy([1 1 0]);
       end
       if alpha_t > alpha_max0
          alpha_t = alpha_max0;
          [phi_psi_alpha_max0,g_max0,~] = feval(f,x_0 + alpha_max0*d,[1 1 0]);
          phi_psi_prim_alpha_max0=g_max0'*d;
          phi_psi_prim_alpha_t = phi_psi_prim_alpha_max0;
          phi_psi_alpha_t = phi_psi_alpha_max0;
          g_t = g_max0;
          evalNumbers.incrementBy([1 1 0]);
       end
       
       if (isBracketed && abs(alpha_l-alpha_u) < workPrec*alpha_max)
           
           outT = alpha_l;
           outX = (x_0+alpha_l*d)';
           outVal = phi_psi_alpha_l;
           outGr = g_l;
           return;
       end
       
       if (alpha_t > 0 && ...
           phi_psi_alpha_t <= f_at_x_0 + mu*alpha_t*phi_prim_0 && ...
           abs(phi_psi_prim_alpha_t) <= eta*abs(phi_prim_0))       
           
           outT = alpha_t;
           outX = (x_0+alpha_t*d)';
           outVal = phi_psi_alpha_t;
           outGr = g_t;
           return;
       end
           
       if (stage1 && phi_psi_alpha_t - f_at_x_0 - mu*alpha_t*phi_prim_0 <= 0 && phi_psi_prim_alpha_t - mu*phi_prim_0 > 0)
          stage1 = false;
       end
      
       if(stage1 && phi_psi_alpha_t <= phi_psi_alpha_l &&  phi_psi_alpha_t > f_at_x_0 + mu*alpha_t*phi_prim_0)
           phi_psi_alpha_t = phi_psi_alpha_t - f_at_x_0 - mu*alpha_t*phi_prim_0;
           phi_psi_alpha_l = phi_psi_alpha_l - f_at_x_0 - mu*alpha_l*phi_prim_0;
           phi_psi_alpha_u = phi_psi_alpha_u - f_at_x_0 - mu*alpha_u*phi_prim_0;
           phi_psi_prim_alpha_t = phi_psi_prim_alpha_t - mu*phi_prim_0;
           phi_psi_prim_alpha_l = phi_psi_prim_alpha_l - mu*phi_prim_0;
           phi_psi_prim_alpha_u = phi_psi_prim_alpha_u - mu*phi_prim_0;
           alpha_t_new = nextSafeguard();           
           %generateNextInterval();
           [alpha_u,alpha_l,phi_psi_alpha_u,phi_psi_alpha_l,phi_psi_prim_alpha_u,phi_psi_prim_alpha_l] = generateNextInterval2(alpha_u,alpha_l,alpha_t,phi_psi_alpha_u,phi_psi_alpha_l,phi_psi_alpha_t,phi_psi_prim_alpha_u,phi_psi_prim_alpha_l,phi_psi_prim_alpha_t);
           phi_psi_alpha_t = phi_psi_alpha_t + f_at_x_0 + mu*alpha_t*phi_prim_0;
           phi_psi_alpha_l = phi_psi_alpha_l - f_at_x_0 + mu*alpha_l*phi_prim_0;
           phi_psi_alpha_u = phi_psi_alpha_u - f_at_x_0 + mu*alpha_u*phi_prim_0;
           phi_psi_prim_alpha_t = phi_psi_prim_alpha_t + mu*phi_prim_0;
           phi_psi_prim_alpha_l = phi_psi_prim_alpha_l + mu*phi_prim_0;
           phi_psi_prim_alpha_u = phi_psi_prim_alpha_u + mu*phi_prim_0;
       else
           alpha_t_new = nextSafeguard();
           %generateNextInterval();
           [alpha_u,alpha_l,phi_psi_alpha_u,phi_psi_alpha_l,phi_psi_prim_alpha_u,phi_psi_prim_alpha_l] = generateNextInterval2(alpha_u,alpha_l,alpha_t,phi_psi_alpha_u,phi_psi_alpha_l,phi_psi_alpha_t,phi_psi_prim_alpha_u,phi_psi_prim_alpha_l,phi_psi_prim_alpha_t);
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
           outGr = g_1;
           return;
       end
       
       [phi_psi_alpha_t,g_t,~] = feval(f,x_0 + alpha_t*d,[1 1 0]);
       phi_psi_prim_alpha_t=g_t'*d;       
       evalNumbers.incrementBy([1 1 0]);
              
    end

    function [alpha_t_new] = nextSafeguard()
        if(phi_psi_alpha_t > phi_psi_alpha_l)
            if(alpha_t == alpha_l)                
                theta = 3.*(phi_psi_alpha_l - phi_psi_alpha_t)/(alpha_t - alpha_l + workPrec) + phi_psi_prim_alpha_l + phi_psi_prim_alpha_t;
            else                
                theta = 3.*(phi_psi_alpha_l - phi_psi_alpha_t)/(alpha_t - alpha_l) + phi_psi_prim_alpha_l + phi_psi_prim_alpha_t;
            end
            s = norm([theta,phi_psi_prim_alpha_l,phi_psi_prim_alpha_t],inf);
            gamma = s*sqrt((theta/s).^2 - (phi_psi_prim_alpha_l/s).*(phi_psi_prim_alpha_t/s));
            if (alpha_t < alpha_l) 
                gamma = -gamma;
            end
            p = (gamma - phi_psi_prim_alpha_l) + theta;
            q = ((gamma - phi_psi_prim_alpha_l) + gamma) + phi_psi_prim_alpha_t;
            r = p/q;
            alpha_c = alpha_l + r*(alpha_t - alpha_l);
            
            alpha_q = alpha_l + ((phi_psi_prim_alpha_l/((phi_psi_alpha_l - phi_psi_alpha_t)/(alpha_t-alpha_l)+phi_psi_prim_alpha_l))/2).*(alpha_t - alpha_l);
            
            if abs(alpha_c-alpha_l) < abs(alpha_q-alpha_l)
                alpha_t_new = alpha_c;
            else
                alpha_t_new = 0.5*(alpha_q+alpha_c);
            end
            isBracketed = true;
            refinePoint = true;                       
        else
            if(alpha_t==alpha_l)
                theta = 3*(phi_psi_alpha_l - phi_psi_alpha_t)/(alpha_t - alpha_l + workPrec) + phi_psi_prim_alpha_l + phi_psi_prim_alpha_t;
            else               
                theta = 3*(phi_psi_alpha_l - phi_psi_alpha_t)/(alpha_t - alpha_l) + phi_psi_prim_alpha_l + phi_psi_prim_alpha_t; 
            end
            s = norm([theta,phi_psi_prim_alpha_l,phi_psi_prim_alpha_t],inf);
            gamma = s*sqrt((theta/s).^2 - (phi_psi_prim_alpha_l/s).*(phi_psi_prim_alpha_t/s));
            if (alpha_t > alpha_l) 
                gamma = -gamma;
            end

            p = (gamma - phi_psi_prim_alpha_t) + theta;
            q = ((gamma - phi_psi_prim_alpha_t) + gamma) + phi_psi_prim_alpha_l;
            r = p/q;
            alpha_c = alpha_t + r*(alpha_l - alpha_t);
            
            alpha_s = alpha_t + (phi_psi_prim_alpha_t/(phi_psi_prim_alpha_t-phi_psi_prim_alpha_l))*(alpha_l - alpha_t);   
            
            if phi_psi_prim_alpha_t*phi_psi_prim_alpha_l < 0
                if abs(alpha_c-alpha_t) > abs(alpha_s-alpha_t)
                   alpha_t_new = alpha_c;
                else
                   alpha_t_new = alpha_s;
                end               
                isBracketed = true;
                refinePoint = false;           
            else
                if abs(phi_psi_prim_alpha_t) < abs(phi_psi_prim_alpha_l)
                    if alpha_t == alpha_l
                        theta = 3*(phi_psi_alpha_l - phi_psi_alpha_t)/(alpha_t - alpha_l+workPrec) + phi_psi_prim_alpha_l + phi_psi_prim_alpha_t;
                    else                       
                        theta = 3*(phi_psi_alpha_l - phi_psi_alpha_t)/(alpha_t - alpha_l) + phi_psi_prim_alpha_l + phi_psi_prim_alpha_t; 
                    end
                    s = norm([theta,phi_psi_prim_alpha_l,phi_psi_prim_alpha_t],inf);
                    gamma = s*sqrt(max(0.,(theta/s)^2 - (phi_psi_prim_alpha_l/s)*(phi_psi_prim_alpha_t/s)));
                    if (alpha_t > alpha_l) 
                        gamma = -gamma;
                    end
                    p = (gamma - phi_psi_prim_alpha_t) + theta;
                    q = (gamma + (phi_psi_prim_alpha_l - phi_psi_prim_alpha_t)) + gamma;
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
                            theta = 3*(phi_psi_alpha_t - phi_psi_alpha_u)/(alpha_u - alpha_t+workPrec) + phi_psi_prim_alpha_u + phi_psi_prim_alpha_t;
                        else
                            theta = 3*(phi_psi_alpha_t - phi_psi_alpha_u)/(alpha_u - alpha_t) + phi_psi_prim_alpha_u + phi_psi_prim_alpha_t;
                        end
                        s = norm([theta,phi_psi_prim_alpha_u,phi_psi_prim_alpha_t],inf);
                        gamma = s*sqrt((theta/s)^2 - (phi_psi_prim_alpha_u/s)*(phi_psi_prim_alpha_t/s));
                        if (alpha_t > alpha_u) 
                            gamma = -gamma;
                        end
                        p = (gamma - phi_psi_prim_alpha_t) + theta;
                        q = ((gamma - phi_psi_prim_alpha_t) + gamma) + phi_psi_prim_alpha_u;
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

    function [] = generateNextInterval()
        
       if (phi_psi_alpha_t > phi_psi_alpha_l)
           alpha_u = alpha_t;
           phi_psi_alpha_u = phi_psi_alpha_t;
           phi_psi_prim_alpha_u = phi_psi_prim_alpha_t;
       else
           v = phi_psi_prim_alpha_t*(alpha_t-alpha_l);
           if (v < 0)
               alpha_l = alpha_t;
               phi_psi_alpha_l = phi_psi_alpha_t;
               phi_psi_prim_alpha_l = phi_psi_prim_alpha_t;
           else
               alpha_u = alpha_l;
               alpha_l = alpha_t;
               phi_psi_alpha_u = phi_psi_alpha_l;
               phi_psi_prim_alpha_u = phi_psi_prim_alpha_l;
               phi_psi_alpha_l = phi_psi_alpha_t;
               phi_psi_prim_alpha_l = phi_psi_prim_alpha_t;
           end
       end
    end
    function [alpha_u_new,alpha_l_new,phi_psi_alpha_u_new,phi_psi_alpha_l_new,phi_psi_prim_alpha_u_new,phi_psi_prim_alpha_l_new] = generateNextInterval2(alpha_u,alpha_l,alpha_t,phi_psi_alpha_u,phi_psi_alpha_l,phi_psi_alpha_t,phi_psi_prim_alpha_u,phi_psi_prim_alpha_l,phi_psi_prim_alpha_t)
       
       alpha_u_new=alpha_u;
       alpha_l_new=alpha_l;
       phi_psi_alpha_u_new=phi_psi_alpha_u;
       phi_psi_alpha_l_new=phi_psi_alpha_l;
       phi_psi_prim_alpha_u_new = phi_psi_prim_alpha_u;
       phi_psi_prim_alpha_l_new = phi_psi_prim_alpha_l;
       
       if (phi_psi_alpha_t > phi_psi_alpha_l)
           alpha_u_new = alpha_t;
           phi_psi_alpha_u_new = phi_psi_alpha_t;
           phi_psi_prim_alpha_u_new = phi_psi_prim_alpha_t;
       else
           v = phi_psi_prim_alpha_t*(alpha_t-alpha_l);
           if (v < 0)
               alpha_l_new = alpha_t;
               phi_psi_alpha_l_new = phi_psi_alpha_t;
               phi_psi_prim_alpha_l_new = phi_psi_prim_alpha_t;
           else
               alpha_u_new = alpha_l;
               alpha_l_new = alpha_t;
               phi_psi_alpha_u_new = phi_psi_alpha_l;
               phi_psi_prim_alpha_u_new = phi_psi_prim_alpha_l;
               phi_psi_alpha_l_new = phi_psi_alpha_t;
               phi_psi_prim_alpha_l_new = phi_psi_prim_alpha_t;
           end
       end
    end
        
end

